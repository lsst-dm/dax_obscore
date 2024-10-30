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
from collections.abc import Iterator
from functools import cache
from typing import Any, cast

import astropy.io.votable
import astropy.table
import felis.datamodel
import pyarrow
import sqlalchemy
import yaml
from lsst.daf.butler import Butler, DataCoordinate, Dimension, Registry, ddl
from lsst.daf.butler.formatters.parquet import arrow_to_numpy
from lsst.daf.butler.registry.obscore import (
    ExposureRegionFactory,
    ObsCoreSchema,
    RecordFactory,
    SpatialObsCorePlugin,
)
from lsst.resources import ResourcePath
from lsst.sphgeom import Region
from lsst.utils.logging import getLogger
from numpy import ma
from pyarrow import RecordBatch, Schema
from pyarrow import Table as ArrowTable
from pyarrow.csv import CSVWriter, WriteOptions
from pyarrow.parquet import ParquetWriter

from . import ExporterConfig
from .config import WhereBind

_LOG = getLogger(__name__)

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


@cache
def _get_obscore_schema() -> felis.datamodel.Schema:
    """Read the ObsCore schema definition."""
    obscore_defn = ResourcePath("resource://lsst.dax.obscore/configs/obscore_nominal.yaml").read()
    obscore_data = yaml.safe_load(obscore_defn)
    schema: felis.datamodel.Schema = felis.datamodel.Schema.model_validate(obscore_data)
    return schema


class _BatchCollector:
    """Helper class to collect records data before making a record batch."""

    def __init__(self, schema: Schema):
        self.schema = schema
        self.batch: list[list] = [[] for column in self.schema.names]
        self.size = 0

    def add_to_batch(self, data: dict[str, Any]) -> None:
        """Add new row to a batch.

        Parameters
        ----------
        data : `dict` [`str`, `~typing.Any`]
            New row.

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
        """Write next buffer to output.

        Parameters
        ----------
        buffer : `bytes`
            Buffer to write.

        Returns
        -------
        size : `int`
            The size of the buffer written.
        """
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
        self._visit_regions: dict[str, dict[int, Region]] = {}
        # Maps instrument+visit+detector to a region
        self._visit_detector_regions: dict[str, dict[tuple[int, int], Region]] = {}
        # Maps instrument and exposure ID to a visit ID
        self._exposure_to_visit: dict[str, dict[int, int]] = {}

    def exposure_region(self, dataId: DataCoordinate) -> Region | None:
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
        if str(detector_dimension) in dataId:
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
            for record_batch, _ in self._make_record_batches(self.config.batch_size):
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
                for record_batch, _ in self._make_record_batches(self.config.batch_size):
                    writer.write_batch(record_batch)

    def to_votable(self, limit: int | None = None) -> astropy.io.votable.tree.VOTableFile:
        """Run the export and return the results as a VOTable instance.

        Parameters
        ----------
        limit : `int` or `None`, optional
            Maximum number of records to return. If `None` there is no limit.

        Returns
        -------
        votable : `astropy.io.votable.tree.VOTableFile`
            The resulting matches as a VOTable.
        """
        # Get the (possibly cached) ObsCore schema.
        schema = _get_obscore_schema()

        tables = schema.tables
        if len(tables) != 1:
            raise RuntimeError("More than one table defined in ObsCore schema")
        obscore_columns = {column.name: column for column in tables[0].columns}

        votable = astropy.io.votable.tree.VOTableFile()
        resource = astropy.io.votable.tree.Resource()
        votable.resources.append(resource)

        fields = []
        for arrow_field in self.schema:
            if arrow_field.name in obscore_columns:
                ffield = obscore_columns[arrow_field.name]
                votable_datatype = felis.datamodel.FelisType.felis_type(ffield.datatype.value).votable_name
                field = astropy.io.votable.tree.Field(
                    votable,
                    name=ffield.name,
                    datatype=votable_datatype,
                    arraysize=ffield.votable_arraysize,
                    unit=ffield.ivoa_unit,
                    ucd=ffield.ivoa_ucd,
                    utype=ffield.votable_utype,
                )
                fields.append(field)
            elif arrow_field.name == "em_filter_name":
                # Non-standard but part of internal standard schema.
                field = astropy.io.votable.tree.Field(
                    votable,
                    name="em_filter_name",
                    datatype="char",
                    arraysize="*",
                )
                fields.append(field)
            else:
                # This must be a non-standard field. Attempt to add it.
                kwargs = {}
                datatype = ""
                if (
                    arrow_field.type.equals(pyarrow.int64())
                    or arrow_field.type.equals(pyarrow.int32())
                    or arrow_field.type.equals(pyarrow.int16())
                ):
                    datatype = "int"
                elif arrow_field.type.equals(pyarrow.string()):
                    datatype = "char"
                    kwargs["arraysize"] = "*"
                elif arrow_field.type.equals(pyarrow.float64()):
                    datatype = "double"
                else:
                    raise RuntimeError(f"Could not handle unrecognized column {arrow_field.name}")
                field = astropy.io.votable.tree.Field(
                    votable,
                    name=arrow_field.name,
                    datatype=datatype,
                    **kwargs,
                )
                fields.append(field)

        table0 = astropy.io.votable.tree.TableElement(votable)
        resource.tables.append(table0)
        table0.fields.extend(fields)

        chunks = []
        n_rows = 0
        overflow = False
        for record_batch, overflow in self._make_record_batches(self.config.batch_size, limit=limit):
            table = ArrowTable.from_batches([record_batch])
            chunk = arrow_to_numpy(table)
            n_rows += len(chunk)
            chunks.append(chunk)

        # Report any overflow.
        query_status = "OVERFLOW" if overflow else "OK"
        info = astropy.io.votable.tree.Info(name="QUERY_STATUS", value=query_status)
        resource.infos.append(info)

        # Combine all the chunks.
        if chunks:
            table0.array = ma.hstack(chunks)

        # Write the output file.
        _LOG.info("Got %d result%s%s", n_rows, "" if n_rows == 1 else "s", " (overflow)" if overflow else "")
        return votable

    def to_votable_file(self, output: str, limit: int | None = None) -> None:
        """Export Butler datasets as ObsCore data model in VOTable format.

        Parameters
        ----------
        output : `str`
            Location of the output file.
        limit : `int` or `None`
            Limit on number of records to return. `None` means no limit.
        """
        votable = self.to_votable(limit=limit)
        votable.to_xml(output)

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

    def _make_record_batches(
        self, batch_size: int = 10_000, limit: int | None = None
    ) -> Iterator[tuple[RecordBatch, bool]]:
        """Generate batches of records to save to a file.

        Yields the batches and a flag indicating whether an overflow condition
        was hit.
        """
        batch = _BatchCollector(self.schema)

        # Set overflow flag.
        overflow = False

        collections: Any = self.config.collections
        if not collections:
            raise ValueError("No collections specified. Querying all collections is not allowed.")

        if limit is not None:
            if limit == 0:
                # Return immediately since no records requested.
                return
            # Always ask for one extra to allow overflow detection.
            limit = abs(limit) + 1

        for dataset_type_name in self.config.dataset_types:
            _LOG.verbose("Querying datasets for dataset type %s [limit=%s]", dataset_type_name, limit)
            where_clauses = self.config.dataset_type_constraints.get(dataset_type_name, [self.config.where])
            if not where_clauses:
                # Want an empty default to match everything.
                where_clauses = [WhereBind(where="")]

            with self.butler.query() as query:
                for where_clause in where_clauses:
                    if where_clause.extra_dims:
                        query = query.join_dimensions(where_clause.extra_dims)

                    if where_clause.where:
                        _LOG.verbose("Processing query with constraint %s", where_clause)
                        query = query.where(where_clause.where, bind=where_clause.bind)

                    refs = query.datasets(dataset_type_name, collections=collections, find_first=True)

                    if limit is not None:
                        refs = refs.limit(limit)

                    # need dimension records
                    count = 0
                    for ref in refs.with_dimension_records():
                        dataId = ref.dataId
                        _LOG.debug("New record, dataId=%s", dataId.mapping)
                        # _LOG.debug("New record, records=%s", dataId.records)

                        record = self.record_factory(ref)
                        if record is None:
                            continue

                        count += 1
                        if limit is not None and count == limit:
                            # Hit the +1 so should not add this to the batch.
                            _LOG.debug("Got one more than requested limit so dropping final record.")
                            overflow = True
                            break

                        batch.add_to_batch(record)
                        if batch.size >= batch_size:
                            _LOG.debug("Saving next record batch, size=%s", batch.size)
                            yield (batch.make_record_batch(), overflow)

                    if limit is not None:
                        limit -= count
                    if overflow:
                        # We counted one too many so adjust for the log
                        # message.
                        count -= 1

                    _LOG.info("Copied %d records from dataset type %s", count, dataset_type_name)

                    if overflow:
                        # No more queries need to run.
                        # This breaks out one level of nesting.
                        break

                if overflow:
                    # Stop further dataset type queries.
                    break

        # Final batch if anything is there
        if batch.size > 0:
            _LOG.debug("Saving final record batch, size=%s", batch.size)
            yield (batch.make_record_batch(), overflow)
