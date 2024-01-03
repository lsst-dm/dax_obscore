# This file is part of dax_obscore.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (http://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import annotations

__all__ = ["ObscoreExporter"]

import contextlib
import io
import logging
from typing import TYPE_CHECKING, Any, Dict, Iterator, List, Optional, Tuple, cast

import pyarrow
import sqlalchemy
from lsst.daf.butler import Butler, DataCoordinate, Dimension, Registry, ddl
from lsst.daf.butler.registry.obscore import (
    ExposureRegionFactory,
    ObsCoreSchema,
    RecordFactory,
    SpatialObsCorePlugin,
)
from lsst.daf.butler.registry.queries import SqlQueryBackend
from lsst.daf.butler.registry.sql_registry import SqlRegistry
from lsst.sphgeom import Region
from pyarrow import RecordBatch, Schema
from pyarrow.csv import CSVWriter, WriteOptions
from pyarrow.parquet import ParquetWriter

from . import ExporterConfig

if TYPE_CHECKING:
    from lsst.daf.butler.registry.queries import SqlQueryContext

_LOG = logging.getLogger(__name__)

# Map few standard Python types to pyarrow types
_PYARROW_TYPE = {
    sqlalchemy.Boolean: pyarrow.bool_(),
    sqlalchemy.SmallInteger: pyarrow.int16(),
    sqlalchemy.Integer: pyarrow.int32(),
    sqlalchemy.BigInteger: pyarrow.int64(),
    sqlalchemy.Float: pyarrow.float64(),
    sqlalchemy.String: pyarrow.string(),
    sqlalchemy.Text: pyarrow.string(),
}


class _BatchCollector:
    """Helper class to collect records data before making a record batch."""

    def __init__(self, schema: Schema):
        self.schema = schema
        self.batch: List[List] = [[] for column in self.schema.names]
        self.size = 0

    def add_to_batch(self, data: Dict[str, Any]) -> None:
        """Add new row to a batch.

        Notes
        -----
        `data` dictionary is updated in place for efficiency.
        """
        for i, column in enumerate(self.schema.names):
            value = data.pop(column, None)
            self.batch[i].append(value)
        self.size += 1

        # watch for unknown columns
        if data:
            columns = set(data.keys())
            raise ValueError(f"Unexpected column names: {columns}")

    def make_record_batch(self) -> RecordBatch:
        """Make a record batch out of accumulated data, and reset."""
        if self.size == 0:
            return None

        # make pyarrow batch out of collected data
        batch = pyarrow.record_batch(self.batch, self.schema)

        # reset to empty
        self.batch = [[] for column in self.schema.names]
        self.size = 0

        return batch


class _CSVFile(io.BufferedWriter):
    r"""File object that intercepts output data and does some editing.

    Parameters
    ----------
    fname : `str`
        Name for the output CSV file.
    null_value : `bytes`
        Value to insert into empty (NULL) cells.
    sep_in : `bytes`
        Field delimiter that appears in input data. Should be selected to be
        something that never appears in actual filed values, e.g. non-printable
        ASCII value.
    sep_out : `bytes`
        Replacement value for field separator, e.g. comma.

    Notes
    -----
    This is a dirty hack to allow writing "\N" to CSV for NULL values instead
    of empty cells. Should be removed when "null_string" can be specified in
    WriteOptions to CSVWriter class.
    """

    def __init__(self, fname: str, null_value: bytes, sep_in: bytes, sep_out: bytes):
        rawfile = open(fname, "wb", buffering=0)
        super().__init__(rawfile)
        self.null_value = null_value
        self.sep_in = sep_in
        self.sep_out = sep_out
        self.buffer: bytes = b""

    def write(self, buffer: bytes) -> int:  # type: ignore
        """Write next buffer to output."""
        self.buffer += buffer
        self._process_buffer()
        return len(buffer)

    def _process_buffer(self, final: bool = False) -> None:
        """Process current buffer contents and write processed data.

        Parameters
        ----------
        final : `bool`
            If True then do not expect any more input data.
        """
        # Replace empty fields with <null_value>.
        while self.buffer:
            line, sep, remainder = self.buffer.partition(b"\n")
            if not sep and not final:
                # no new line, wait for more input
                break
            # keep everything after new line for next round
            self.buffer = remainder
            # Replace empty cells, and separators
            line = self.sep_out.join(cell if cell else self.null_value for cell in line.split(self.sep_in))
            super().write(line + sep)

    def close(self) -> None:
        """Finalize writing."""
        self._process_buffer(final=True)
        super().close()


class _ExposureRegionFactory(ExposureRegionFactory):
    """Find exposure region from a matching visit dimensions records."""

    def __init__(self, registry: Registry):
        self.registry = registry
        self.universe = registry.dimensions

        # Maps instrument and visit ID to a region
        self._visit_regions: Dict[str, Dict[int, Region]] = {}
        # Maps instrument+visit+detector to a region
        self._visit_detector_regions: Dict[str, Dict[Tuple[int, int], Region]] = {}
        # Maps instrument and exposure ID to a visit ID
        self._exposure_to_visit: Dict[str, Dict[int, int]] = {}

    def exposure_region(self, dataId: DataCoordinate, context: SqlQueryContext) -> Optional[Region]:
        # Docstring is inherited from a base class.
        registry = self.registry
        instrument = cast(str, dataId["instrument"])

        exposure_to_visit = self._exposure_to_visit.get(instrument)
        if exposure_to_visit is None:
            self._exposure_to_visit[instrument] = exposure_to_visit = {}
            # Read complete relation between visits and exposures. There could
            # be multiple visits defined per exposure, but they are supposed to
            # have the same region, so we take one of them at random.
            records = registry.queryDimensionRecords("visit_definition", instrument=instrument)
            for record in records:
                exposure_to_visit[record.exposure] = record.visit
            _LOG.debug("read %d exposure-to-visit records", len(exposure_to_visit))

        # map exposure to a visit
        exposure = cast(int, dataId["exposure"])
        visit = exposure_to_visit.get(exposure)
        if visit is None:
            return None

        universe = self.universe
        detector_dimension = cast(Dimension, universe["detector"])
        if detector_dimension in dataId:
            visit_detector_regions = self._visit_detector_regions.get(instrument)

            if visit_detector_regions is None:
                self._visit_detector_regions[instrument] = visit_detector_regions = {}

                # Read all visits, there is a chance we need most of them
                # anyways, and trying to filter by dataset type and collection
                # makes it much slower.
                records = registry.queryDimensionRecords("visit_detector_region", instrument=instrument)
                for record in records:
                    visit_detector_regions[(record.visit, record.detector)] = record.region
                _LOG.debug("read %d visit-detector regions", len(visit_detector_regions))

            detector = cast(int, dataId["detector"])
            return visit_detector_regions.get((visit, detector))

        else:
            visit_regions = self._visit_regions.get(instrument)

            if visit_regions is None:
                self._visit_regions[instrument] = visit_regions = {}

                # Read all visits, there is a chance we need most of them
                # anyways, and trying to filter by dataset type and collection
                # makes it much slower.
                records = registry.queryDimensionRecords("visit", instrument=instrument)
                for record in records:
                    visit_regions[record.id] = record.region
                _LOG.debug("read %d visit regions", len(visit_regions))

            return visit_regions.get(visit)


class ObscoreExporter:
    """Class for extracting and exporting of the datasets in ObsCore format.

    Parameters
    ----------
    butler : `lsst.daf.butler.Butler`
        Data butler.
    config : `lsst.dax.obscore.ExporterConfig`
        Exporter configuration.
    """

    def __init__(self, butler: Butler, config: ExporterConfig):
        self.butler = butler
        self.config = config

        # Build schema for a table from ObsCoreSchema and plugins.
        spatial_plugins = SpatialObsCorePlugin.load_plugins(config.spatial_plugins, None)
        schema = ObsCoreSchema(config=config, spatial_plugins=spatial_plugins)

        self.schema = self._make_schema(schema.table_spec)

        exposure_region_factory = _ExposureRegionFactory(self.butler.registry)
        universe = self.butler.dimensions
        self.record_factory = RecordFactory(
            config, schema, universe, spatial_plugins, exposure_region_factory
        )

    def to_parquet(self, output: str) -> None:
        """Export Butler datasets as ObsCore Data Model in parquet format.

        Parameters
        ----------
        output : `str`
            Location of the output file.
        """
        compression = self.config.parquet_compression
        with ParquetWriter(output, self.schema, compression=compression) as writer:
            for record_batch in self._make_record_batches(self.config.batch_size):
                writer.write_batch(record_batch)

    def to_csv(self, output: str) -> None:
        """Export Butler datasets as ObsCore Data Model in CSV format.

        Parameters
        ----------
        output : `str`
            Location of the output file.
        """
        # Use Unit Separator (US) = 0x1F as field delimiter to avoid potential
        # issue with commas appearing in actual values.
        options = WriteOptions(delimiter="\x1f")
        null_string = self.config.csv_null_string.encode()
        with contextlib.closing(_CSVFile(output, null_string, sep_in=b"\x1f", sep_out=b",")) as file:
            with CSVWriter(file, self.schema, write_options=options) as writer:
                for record_batch in self._make_record_batches(self.config.batch_size):
                    writer.write_batch(record_batch)

    def _make_schema(self, table_spec: ddl.TableSpec) -> Schema:
        """Create schema definition for output data.

        Parameters
        ----------
        table_spec : `ddl.TableSpec`

        Returns
        -------
        schema : `pyarrow.Schema`
            Schema definition.
        """
        schema: list[tuple] = []
        for field_spec in table_spec.fields:
            try:
                pyarrow_type = _PYARROW_TYPE[field_spec.dtype]
            except KeyError:
                raise TypeError(
                    f"Unexpected type of column column={field_spec.name}, value={field_spec.dtype}"
                ) from None
            schema.append((field_spec.name, pyarrow_type))

        return pyarrow.schema(schema)

    def _make_record_batches(self, batch_size: int = 10_000) -> Iterator[RecordBatch]:
        """Generate batches of records to save to a file."""
        batch = _BatchCollector(self.schema)

        collections: Any = self.config.collections
        if not collections:
            collections = ...

        # Have to use non-public Registry interface.
        registry = self.butler._registry
        assert isinstance(registry, SqlRegistry), "Registry must be SqlRegistry"
        backend = SqlQueryBackend(registry._db, registry._managers)

        context = backend.context()
        for dataset_type_name in self.config.dataset_types:
            _LOG.debug("Reading data for dataset %s", dataset_type_name)
            refs = registry.queryDatasets(dataset_type_name, collections=collections, where=self.config.where)

            # need dimension records
            refs = refs.expanded()
            count = 0
            for ref in refs:
                dataId = ref.dataId
                _LOG.debug("New record, dataId=%s", dataId.full)
                # _LOG.debug("New record, records=%s", dataId.records)

                record = self.record_factory(ref, context)
                if record is None:
                    continue

                count += 1

                batch.add_to_batch(record)
                if batch.size >= batch_size:
                    _LOG.debug("Saving next record batch, size=%s", batch.size)
                    yield batch.make_record_batch()

            _LOG.info("Copied %d records from dataset type %s", count, dataset_type_name)

        # Final batch if anything is there
        if batch.size > 0:
            _LOG.debug("Saving final record batch, size=%s", batch.size)
            yield batch.make_record_batch()
