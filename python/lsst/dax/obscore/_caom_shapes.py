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

__all__ = ["require_caom2"]

try:
    import caom2
except ImportError:
    caom2 = None

_CAOM2_IMPORT_MESSAGE = (
    "The caom2 package is required for CAOM export but is not installed. "
    "Install it with 'pip install lsst-dax-obscore[caom]'."
)


def require_caom2() -> None:
    """Check that the optional ``caom2`` dependency is importable.

    Raises
    ------
    ImportError
        Raised if ``caom2`` is not installed.
    """
    if caom2 is None:
        raise ImportError(_CAOM2_IMPORT_MESSAGE)
