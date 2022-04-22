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

import logging
from typing import Any, Iterator, List, Union

import pyarrow
from lsst.daf.butler import Butler
from pyarrow import parquet, RecordBatch, Schema

from . import DatasetTypeConfig, ExporterConfig


_LOG = logging.getLogger(__name__)

# List of standard columns in output file. Extra columns can be added via
# `extra_columns` parameters in configuration. This should include at least
# all mandatory columns defined in ObsCore note (revision 1.1, Appendix B)
_STATIC_SCHEMA = (
    ("dataproduct_type", pyarrow.string()),
    ("dataproduct_subtype", pyarrow.string()),
    ("calib_level", pyarrow.int8()),
    ("target_name", pyarrow.string()),
    ("obs_id", pyarrow.string()),
    ("obs_collection", pyarrow.string()),
    ("obs_publisher_did", pyarrow.string()),
    ("access_url", pyarrow.string()),
    ("access_format", pyarrow.string()),
    ("s_ra", pyarrow.float64()),
    ("s_dec", pyarrow.float64()),
    ("s_fow", pyarrow.float64()),
    ("s_region", pyarrow.string()),
    ("s_resolution", pyarrow.float64()),
    ("s_xel1", pyarrow.int16()),
    ("s_xel2", pyarrow.int16()),
    ("t_xel", pyarrow.int16()),
    ("t_min", pyarrow.float64()),
    ("t_max", pyarrow.float64()),
    ("t_exptime", pyarrow.float64()),
    ("t_resolution", pyarrow.float64()),
    ("em_xel", pyarrow.int16()),
    ("em_min", pyarrow.float64()),
    ("em_max", pyarrow.float64()),
    ("em_res_power", pyarrow.float64()),
    ("o_ucd", pyarrow.string()),
    ("pol_xel", pyarrow.int16()),
    ("instrument_name", pyarrow.string()),
)


_PYARROW_TYPE = {
    bool: pyarrow.bool_(),
    int: pyarrow.int64(),
    float: pyarrow.float64(),
    str: pyarrow.string(),
}


class ObscoreExporter:
    """Class for extracting and exporting of the datasets in ObsCore format.

    Parameters
    ----------
    butler : `lsst.daf.butler.Butler`
        Data butler.
    config : `lsst.daf.butler.Config`
    """

    def __init__(self, butler: Butler, config: ExporterConfig):
        self.butler = butler
        self.config = config
        self.schema = self._make_schema(config)

    def to_parquet(self, output: str) -> None:
        """Export Butler datasets as ObsCore Data Model in parquet format.

        Parameters
        ----------
        output : `str`
            Location of the output file.
        """

        with parquet.ParquetWriter(output, self.schema) as writer:
            for record_batch in self._make_record_batches(
                self.config.batch_size
            ):
                writer.write_batch(record_batch)

    def _make_schema(self, config: ExporterConfig) -> Schema:
        """Create schema definition for output data.

        Returns
        -------
        schema : `pyarrow.Schema`
            Schema definition.
        """
        schema = list(_STATIC_SCHEMA)

        columns = set(col[0] for col in schema)

        all_configs: List[Union[ExporterConfig, DatasetTypeConfig]] = [config]
        if config.dataset_types:
            all_configs += config.dataset_types
        for cfg in all_configs:
            if cfg.extra_columns:
                for col_name, col_value in cfg.extra_columns.items():
                    if col_name in columns:
                        continue
                    col_type = _PYARROW_TYPE.get(type(col_value))
                    if col_type is None:
                        raise TypeError(
                            f"Unexpected type in extra_columns: column={col_name}, value={col_value:r}"
                        )
                    schema.append((col_name, col_type))
                    columns.add(col_name)

        return pyarrow.schema(schema)

    def _make_record_batches(self, batch_size: int = 10_000) -> Iterator[RecordBatch]:
        """Generate batches of records to save to a file."""
        collections: Any = self.config.collections
        if not collections:
            collections = ...

        registry = self.butler.registry
        for dataset_config in self.config.dataset_types:
            result = registry.queryDatasets(
                dataset_config.name, collections=collections
            )
            _LOG.debug(
                "Dataset type %s returned %s datasets",
                dataset_config.name,
                result.count(),
            )
