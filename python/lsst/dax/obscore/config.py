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

__all__ = ["ExporterConfig", "WhereBind"]

from typing import Any

from lsst.daf.butler.registry.obscore import ObsCoreConfig
from pydantic import BaseModel, ConfigDict, Field


class WhereBind(BaseModel):
    """A where expression with associated bind parameters."""

    model_config = ConfigDict(frozen=True)

    where: str = ""
    """User expression to restrict the output."""
    bind: dict[str, Any] = Field(default_factory=dict)
    """Bind values specified in the ``where`` expression."""

    @classmethod
    def combine(cls, wheres: list[WhereBind]) -> WhereBind:
        """Combine multiple clauses into a single where expression.

        Parameters
        ----------
        wheres : `list` [ `WhereBind`]
            The user expressions to combine.
        """
        where = " AND ".join(w.where for w in wheres)
        bind = {}
        for w in wheres:
            bind.update(w.bind)
        return cls(where=where, bind=bind)


class ExporterConfig(ObsCoreConfig):
    """Complete configuration for ObscoreExporter."""

    where: WhereBind = Field(default_factory=WhereBind)
    """Default user expression to restrict the output. Not used if
    per dataset type user expression is provided."""

    dataset_type_constraints: dict[str, list[WhereBind]] = Field(default_factory=dict)
    """Specific user expressions for a given dataset type. If a dataset type
    is not specified here the default ``where`` will be used."""

    batch_size: int = 10_000
    """Number of records in a pyarrow RecordBatch"""

    parquet_compression: str = "snappy"
    """Compression method for parquet files"""

    csv_null_string: str = r"\N"
    """Value to use for NULLs in CSV output."""
