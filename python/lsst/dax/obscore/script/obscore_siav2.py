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

__all__ = ["obscore_siav2"]


def obscore_siav2(
    repo: str,
    config: str,
    pos: str | None,
) -> None:
    """Query butler using SIAv2 parameters and return VOTable.

    Parameters
    ----------
    repo : `str`
        URI to the butler repository.
    config : `str`
        Location of the ObsCore configuration file.
    pos : `str` | None
        POS string specifying a region on the sky.
    """
    pass
