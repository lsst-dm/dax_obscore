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

__all__ = ["get_siav2_handler"]

from importlib.metadata import EntryPoint, entry_points
from typing import TYPE_CHECKING

from lsst.utils.introspection import get_full_type_name

if TYPE_CHECKING:
    from .siav2 import SIAv2Handler


def _get_siav2_entry_points() -> dict[str, EntryPoint]:
    plugins = entry_points(group="dax_obscore.siav2")
    return {p.name: p for p in plugins}


def get_siav2_handler(namespace: str) -> type[SIAv2Handler]:
    """Select the correct handler for this universe namespace.

    Parameters
    ----------
    namespace : `str`
        The butler dimension universe namespace.

    Returns
    -------
    handler : `type` [ `lsst.dax.obscore.siav2.SIAv2Handler` ]
        The class of handler suitable for this namespace.
    """
    plugins = _get_siav2_entry_points()
    entry_point = plugins.get(namespace)
    if entry_point is None:
        known = ", ".join(plugins.keys())
        raise RuntimeError(
            f"Unable to find suitable SIAv2 Handler for namespace {namespace} [do understand: {known}]"
        )
    func = entry_point.load()
    handler_type = func()

    # Check that we have the right type. This has to be a deferred load but
    # the code should already be loaded at this point.
    from .siav2 import SIAv2Handler

    if not issubclass(handler_type, SIAv2Handler):
        raise TypeError(
            f"Entry point for universe {namespace} did not return SIAv2Handler. "
            f"Returned {get_full_type_name(handler_type)}"
        )

    return func()
