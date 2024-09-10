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

import math
import unittest

import lsst.sphgeom
from astropy.time import Time
from lsst.daf.butler import Timespan
from lsst.dax.obscore.siav2 import Interval, SIAv2Parameters
from lsst.utils.iteration import ensure_iterable


class SIAv2IntervalTestCase(unittest.TestCase):
    """Test the Interval class."""

    def test_interval(self):
        i1 = Interval(start=3, end=5)
        i2 = Interval(start=4, end=6)
        self.assertTrue(i1.overlaps(i2))
        self.assertEqual(tuple(i1), (3.0, 5.0))

        with self.assertRaises(ValueError):
            Interval(start=3.14, end=-math.inf)

    def test_parsing(self):
        for s, expected in (
            ("-Inf 42", (-math.inf, 42.0)),
            ("-Inf +Inf", (-math.inf, math.inf)),
            ("-5 +Inf", (-5.0, math.inf)),
            ("-10 -5", (-10.0, -5.0)),
        ):
            interval = Interval.from_string(s)
            self.assertEqual(tuple(interval), expected)

        with self.assertRaises(ValueError):
            Interval.from_string("")
        with self.assertRaises(ValueError):
            Interval.from_string("1 2 3")
        with self.assertRaises(ValueError):
            Interval.from_string("HSC 5")


class SIAv2ParametersTestCase(unittest.TestCase):
    """Test parameter parsing."""

    def test_strings(self):
        # String based items without validation.
        for param in ("pol", "instrument", "collection", "facility", "target", "id"):
            for test in ("HSC", ["LSSTCam", "LATISS"], []):
                kwargs = {param: test}
                p = SIAv2Parameters.from_siav2(**kwargs)
                self.assertCountEqual(getattr(p, param), tuple(ensure_iterable(test)))

            with self.assertRaises(ValueError):
                kwargs = {param: 42}
                SIAv2Parameters.from_siav2(**kwargs)

    def test_pos(self):
        for pos in ("CIRCLE 0 0 10", ["RANGE 1 2 3 4", "POLYGON 12.0 34.0 14.0 35.0 14. 36.0 12.0 35.0"], []):
            p = SIAv2Parameters.from_siav2(pos=pos)
            self.assertEqual(len(p.pos), len(list(ensure_iterable(pos))))
            for r in p.pos:
                # Should already be checkd by Pydantic.
                self.assertIsInstance(r, lsst.sphgeom.Region)

    def test_time(self):
        mjds = [58000.0, 58100.0, 59000.0]
        times = [Time(mjd, format="mjd", scale="utc") for mjd in mjds]

        p = SIAv2Parameters.from_siav2(
            time=[
                f"{mjds[0]} {mjds[1]}",
                f"{mjds[0]} {mjds[2]}",
                f"{mjds[1]}",
                f"-Inf {mjds[1]}",
                f"{mjds[0]} +Inf",
            ]
        )
        parsed = p.time
        for index, expected in (
            (0, (times[0], times[1])),
            (1, (times[0], times[2])),
            (3, (None, times[1])),
            (4, (times[0], None)),
        ):
            self.assertEqual(parsed[index], Timespan(expected[0], expected[1]))
        self.assertEqual(parsed[2], times[1])

        with self.assertRaises(ValueError):
            SIAv2Parameters.from_siav2(time="wrong")

    def test_intervals(self):
        p = SIAv2Parameters.from_siav2(
            band="-Inf 42",
            exptime="0 30",
            fov=["1.0 2.0"],
            spatres="0.1 0.2",
            specrp="10000 20000",
            timeres=["1.0 2.0", "-Inf 1.0"],
        )
        self.assertEqual(p.band[0], Interval(start=-math.inf, end=42.0))
        self.assertEqual(p.exptime[0], Interval(start=0, end=30.0))
        self.assertEqual(p.fov[0], Interval(start=1.0, end=2.0))
        self.assertEqual(p.spatres[0], Interval(start=0.1, end=0.2))
        self.assertEqual(p.specrp[0], Interval(start=10000, end=20000))
        self.assertEqual(p.timeres, (Interval(start=1.0, end=2.0), Interval(start=-math.inf, end=1.0)))

    def test_calib(self):
        p = SIAv2Parameters.from_siav2(calib=(1, 2, 3))
        self.assertCountEqual(p.calib, (1, 2, 3))
        p = SIAv2Parameters.from_siav2(calib=0)
        self.assertCountEqual(p.calib, [0])

        with self.assertRaises(ValueError):
            SIAv2Parameters.from_siav2(calib=4)

    def test_dptype(self):
        p = SIAv2Parameters.from_siav2(dptype=["image", "cube"])
        self.assertEqual(p.dptype, frozenset(["image", "cube"]))
        with self.assertRaises(ValueError):
            SIAv2Parameters.from_siav2(dptype="spectrum")


if __name__ == "__main__":
    unittest.main()