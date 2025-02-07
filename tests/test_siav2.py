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
import os
import unittest

from astropy.time import Time

import lsst.sphgeom
from lsst.daf.butler import Timespan
from lsst.daf.butler.tests.utils import makeTestTempDir, removeTestTempDir
from lsst.dax.obscore.plugins import get_siav2_handler
from lsst.dax.obscore.siav2 import Interval, SIAv2Handler, SIAv2Parameters, siav2_query_from_raw
from lsst.dax.obscore.tests import DaxObsCoreTestMixin
from lsst.utils.iteration import ensure_iterable

TESTDIR = os.path.abspath(os.path.dirname(__file__))


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

    def test_maxrec(self):
        p = SIAv2Parameters.from_siav2(maxrec=5)
        self.assertEqual(p.maxrec, 5)
        with self.assertRaises(ValueError):
            SIAv2Parameters.from_siav2(maxrec=-5)


class SIAv2TestCase(unittest.TestCase, DaxObsCoreTestMixin):
    """Tests of SIAv2 queries."""

    def setUp(self):
        self.root = makeTestTempDir(TESTDIR)
        self.butler = self.make_butler()
        self.butler.import_(filename=os.path.join(TESTDIR, "data", "hsc_gen3.yaml"), without_datastore=True)
        self.config = self.make_export_config()

    def tearDown(self):
        removeTestTempDir(self.root)

    def assertVOTable(self, votable, n_rows) -> None:
        """Check a VOTable."""
        tables = list(votable.iter_tables())
        self.assertEqual(len(tables), 1)
        table0 = tables[0].array
        self.assertEqual(len(table0), n_rows, str(table0))

    def test_query(self):
        """Test that an SIAv2 query completes."""
        config = self.config.model_copy()
        config.batch_size = 3
        for kwargs, expected in (
            ({"instrument": "HSC"}, 68),
            ({"instrument": ("HSC", "LATISS")}, 68),
            ({"instrument": "LATISS"}, 2),  # coadds do not know their instrument.
            ({"pos": "CIRCLE 320.851 -0.3 0.001"}, 28),
            ({"band": "700e-9"}, 33),
            ({"band": ("700e-9", "500e-9 600e-9")}, 68),
            ({"band": ("700e-9", "500e-9 600e-9"), "calib": {2}}, 33),
            ({"time": "56460.0 56460.57"}, 36),
            ({"time": ("56460.0 56460.57", "60000 +Inf")}, 36),
            ({"time": "56460.0 56460.57", "pos": "CIRCLE 320.851 -0.3 0.001"}, 16),
            ({"time": "56460.0 56460.57", "pos": "CIRCLE 320.851 -0.3 0.001", "calib": {2}}, 7),
            ({"maxrec": 5}, 5),
            ({"maxrec": 0}, 0),
            ({"exptime": "0 35"}, 66),
            ({"exptime": ("0 5", "40 +Inf")}, 0),
            ({"calib": {2}}, 33),
            ({"calib": {3}}, 2),
            ({"calib": {1}}, 33),
        ):
            with self.subTest(kwargs=kwargs, expected=expected):
                votable = siav2_query_from_raw(
                    self.butler, self.config, collections=["HSC/runs/ci_hsc", "HSC/raw/all"], **kwargs
                )
                self.assertVOTable(votable, expected)

                # Again with a very small batch size to check that looping
                # works.
                votable = siav2_query_from_raw(
                    self.butler, config, collections=["HSC/runs/ci_hsc", "HSC/raw/all"], **kwargs
                )
                self.assertVOTable(votable, expected)

    def test_entry_point(self):
        """Test that a handler can be returned by namespace."""
        handler = get_siav2_handler("daf_butler")
        self.assertTrue(issubclass(handler, SIAv2Handler))

        with self.assertRaises(RuntimeError):
            get_siav2_handler("unknown")


if __name__ == "__main__":
    unittest.main()
