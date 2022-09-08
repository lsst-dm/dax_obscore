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

from typing import Optional

from lsst.daf.butler.registry.obscore import ObsCoreConfig


class ExporterConfig(ObsCoreConfig):
    """Complete configuration for ObscoreExporter."""

    where: Optional[str] = None
    """User expression to restrict the output. This value can be overridden
    with command line options.
    """

    batch_size: int = 10_000
    """Number of records in a pyarrow RecordBatch"""

    parquet_compression: str = "snappy"
    """Compression method for parquet files"""

    csv_null_string: str = r"\N"
    """Value to use for NULLs in CSV output."""
