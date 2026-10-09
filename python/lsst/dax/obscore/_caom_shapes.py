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

__all__ = ["arcsec_to_degrees", "metres", "region_to_caom_shape", "require_caom2"]

from typing import Any

import astropy.units as u

from lsst.sphgeom import Box, Circle, ConvexPolygon, LonLat, Region
from lsst.utils.logging import getLogger

try:
    import caom2
except ImportError:
    caom2 = None

_LOG = getLogger(__name__)

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


def arcsec_to_degrees(value: float | None) -> float | None:
    """Convert an angle in arcseconds to degrees.

    Parameters
    ----------
    value : `float` or `None`
        Angle in arcseconds. `None` passes through unchanged.

    Returns
    -------
    degrees : `float` or `None`
        The angle in degrees.
    """
    if value is None:
        return None
    return float((value * u.arcsec).to_value(u.deg))


def metres(quantity: u.Quantity) -> float:
    """Convert a length quantity to a plain number of metres.

    Parameters
    ----------
    quantity : `astropy.units.Quantity`
        A quantity convertible to metres.

    Returns
    -------
    value : `float`
        The value in metres.
    """
    return float(quantity.to_value(u.m))


def region_to_caom_shape(region: Region | None) -> Any | None:
    """Convert a sphgeom region to the equivalent CAOM shape.

    Parameters
    ----------
    region : `~lsst.sphgeom.Region` or `None`
        Region to convert. `None` passes through unchanged.

    Returns
    -------
    shape : `caom2.shape.Polygon`, `caom2.shape.Circle`, \
            `caom2.shape.Box` or `None`
        The equivalent CAOM shape, or `None` if the region is `None` or of
        an unsupported type.

    Notes
    -----
    An unsupported region type, such as a union region, issues a warning
    and returns `None` rather than raising, so that a single awkward
    dataset does not abort an entire export.
    """
    if region is None:
        return None
    require_caom2()

    if isinstance(region, ConvexPolygon):
        points = [
            caom2.shape.Point(coord.getLon().asDegrees(), coord.getLat().asDegrees())
            for coord in (LonLat(vector) for vector in region.getVertices())
        ]
        # A caom2 Polygon carries both a simple point list and a sampled
        # MultiPolygon; the latter is explicitly closed, so the first
        # vertex is repeated as the final CLOSE segment.
        vertices = [
            caom2.shape.Vertex(point.cval1, point.cval2, caom2.shape.SegmentType.LINE) for point in points
        ]
        vertices[0].type = caom2.shape.SegmentType.MOVE
        vertices.append(caom2.shape.Vertex(points[0].cval1, points[0].cval2, caom2.shape.SegmentType.CLOSE))
        return caom2.shape.Polygon(points=points, samples=caom2.shape.MultiPolygon(vertices=vertices))

    if isinstance(region, Circle):
        center = LonLat(region.getCenter())
        return caom2.shape.Circle(
            center=caom2.shape.Point(center.getLon().asDegrees(), center.getLat().asDegrees()),
            radius=region.getOpeningAngle().asDegrees(),
        )

    if isinstance(region, Box):
        # Box.getCenter() already returns a LonLat, unlike Circle.
        center = region.getCenter()
        lon = region.getLon()
        lat = region.getLat()
        return caom2.shape.Box(
            center=caom2.shape.Point(center.getLon().asDegrees(), center.getLat().asDegrees()),
            width=(lon.getB() - lon.getA()).asDegrees(),
            height=(lat.getB() - lat.getA()).asDegrees(),
        )

    _LOG.warning("Cannot convert region of type %s to a CAOM shape; omitting position.", type(region))
    return None
