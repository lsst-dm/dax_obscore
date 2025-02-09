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

__all__ = ["ExporterConfig"]

from collections.abc import Iterable
from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field

from lsst.daf.butler.registry.obscore import ObsCoreConfig


class WhereBind(BaseModel):
    """A where expression with associated bind parameters."""

    model_config = ConfigDict(frozen=True)

    where: str = ""
    """User expression to restrict the output."""
    bind: dict[str, Any] = Field(default_factory=dict)
    """Bind values specified in the ``where`` expression."""
    extra_dims: frozenset[str] = Field(default_factory=frozenset)
    """Extra dimensions required to be included in query."""

    @classmethod
    def combine(cls, wheres: list[WhereBind], mode: Literal["AND"] | Literal["OR"] = "AND") -> WhereBind:
        """Combine multiple clauses into a single where expression.

        Parameters
        ----------
        wheres : `list` [ `WhereBind`]
            The user expressions to combine.
        mode : `str`
            Combination mode. Can be ``AND`` or ``OR``.

        Returns
        -------
        combo : `WhereBind`
            A new `WhereBind` representing all the information of the input
            clauses.
        """
        n_wheres = len(wheres)
        if n_wheres == 0:
            raise ValueError("Must provide at least one WhereBind for combination.")
        if n_wheres == 1:
            return wheres[0]
        where = f" {mode} ".join(f"({w.where})" for w in wheres if w.where)
        bind: dict[str, Any] = {}
        extras: set[str] = set()
        for w in wheres:
            extras.update(w.extra_dims)
            # For an empty where clause the bind values are irrelevant.
            if not w.where:
                continue
            # Warn if we are overwriting bind keys.
            duplicates = bind.keys() & w.bind.keys()
            if duplicates:
                # Allow duplicates if they are identical.
                differing_keys = set()
                for dupe in duplicates:
                    if bind[dupe] != w.bind[dupe]:
                        differing_keys.add(dupe)
                if differing_keys:
                    raise ValueError(
                        f"Combining multiple WHERE clauses with reused bind parameters of {differing_keys}."
                    )
            bind.update(w.bind)

        return cls(where=where, bind=bind, extra_dims=extras)


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

    def select_dataset_types(self, dataset_types: Iterable[str]) -> None:
        """Update the configuration to include only these dataset types.

        Parameters
        ----------
        dataset_types : `~collections.abc.Iterable` [ `str` ]
            Names of dataset types to select.
        """
        dataset_type_set = set(dataset_types)
        # Check that configuration has all requested dataset types.
        if not dataset_type_set.issubset(self.dataset_types):
            extras = dataset_type_set - set(self.dataset_types)
            raise ValueError(f"Dataset types {extras} are not defined in configuration file.")
        # Remove dataset types that are not needed.
        self.dataset_types = {
            key: value for key, value in self.dataset_types.items() if key in dataset_type_set
        }
