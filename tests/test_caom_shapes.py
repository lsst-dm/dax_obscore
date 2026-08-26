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

import unittest

import astropy.units as u

import lsst.sphgeom
from lsst.dax.obscore import _caom_shapes


class ImportGuardTestCase(unittest.TestCase):
    """Tests of the optional caom2 import guard."""

    def test_require_caom2_succeeds_when_available(self):
        """require_caom2 is a no-op when the library is installed."""
        self.assertIsNotNone(_caom_shapes.caom2)
        _caom_shapes.require_caom2()

    def test_require_caom2_message_names_the_extra(self):
        """The failure message tells the user how to fix the problem."""
        message = _caom_shapes._CAOM2_IMPORT_MESSAGE
        self.assertIn("lsst-dax-obscore[caom]", message)


class ConversionTestCase(unittest.TestCase):
    """Tests of region and unit conversion."""

    def test_convex_polygon(self):
        """A ConvexPolygon becomes a caom2 Polygon with matching vertices."""
        vertices = [
            lsst.sphgeom.UnitVector3d(lsst.sphgeom.LonLat.fromDegrees(lon, lat))
            for lon, lat in ((10.0, 20.0), (11.0, 20.0), (11.0, 21.0), (10.0, 21.0))
        ]
        region = lsst.sphgeom.ConvexPolygon(vertices)

        shape = _caom_shapes.region_to_caom_shape(region)

        self.assertEqual(len(shape.points), 4)
        self.assertEqual(len(shape.samples.vertices), 5)
        expected = {
            (round(v.getLon().asDegrees(), 6), round(v.getLat().asDegrees(), 6))
            for v in (lsst.sphgeom.LonLat(vec) for vec in region.getVertices())
        }
        actual = {(round(p.cval1, 6), round(p.cval2, 6)) for p in shape.points}
        self.assertEqual(actual, expected)

    def test_circle(self):
        """A sphgeom Circle becomes a caom2 Circle in degrees."""
        center = lsst.sphgeom.UnitVector3d(lsst.sphgeom.LonLat.fromDegrees(30.0, -10.0))
        region = lsst.sphgeom.Circle(center, lsst.sphgeom.Angle.fromDegrees(2.5))

        shape = _caom_shapes.region_to_caom_shape(region)

        self.assertAlmostEqual(shape.center.cval1, 30.0, places=6)
        self.assertAlmostEqual(shape.center.cval2, -10.0, places=6)
        self.assertAlmostEqual(shape.radius, 2.5, places=6)

    def test_box(self):
        """A sphgeom Box becomes a caom2 Box in degrees."""
        region = lsst.sphgeom.Box.fromDegrees(10.0, 20.0, 12.0, 23.0)

        shape = _caom_shapes.region_to_caom_shape(region)

        self.assertAlmostEqual(shape.width, 2.0, places=6)
        self.assertAlmostEqual(shape.height, 3.0, places=6)

    def test_none_region(self):
        """A missing region converts to None rather than raising."""
        self.assertIsNone(_caom_shapes.region_to_caom_shape(None))

    def test_unsupported_region_warns(self):
        """An unsupported region logs a warning and yields None."""
        a = lsst.sphgeom.Circle(
            lsst.sphgeom.UnitVector3d(lsst.sphgeom.LonLat.fromDegrees(0.0, 0.0)),
            lsst.sphgeom.Angle.fromDegrees(1.0),
        )
        b = lsst.sphgeom.Circle(
            lsst.sphgeom.UnitVector3d(lsst.sphgeom.LonLat.fromDegrees(90.0, 0.0)),
            lsst.sphgeom.Angle.fromDegrees(1.0),
        )
        region = lsst.sphgeom.UnionRegion(a, b)

        with self.assertLogs("lsst.dax.obscore._caom_shapes", level="WARNING"):
            self.assertIsNone(_caom_shapes.region_to_caom_shape(region))

    def test_arcsec_to_degrees(self):
        """Pixel scale converts from arcsec to degrees via astropy."""
        self.assertAlmostEqual(
            _caom_shapes.arcsec_to_degrees(0.2), (0.2 * u.arcsec).to_value(u.deg), places=12
        )
        self.assertIsNone(_caom_shapes.arcsec_to_degrees(None))

    def test_metres(self):
        """Length quantities convert to plain metres via astropy."""
        self.assertAlmostEqual(_caom_shapes.metres(2.0 * u.km), 2000.0, places=9)


if __name__ == "__main__":
    unittest.main()
