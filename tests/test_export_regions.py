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

import os
import unittest

import astropy.io.fits as fits
from click.testing import CliRunner

from lsst.daf.butler import DatasetType
from lsst.daf.butler.tests.utils import makeTestTempDir, removeTestTempDir
from lsst.dax.obscore.cli.cmd.commands import obscore as obscore_cli
from lsst.dax.obscore.script import obscore_export_regions
from lsst.dax.obscore.tests import DaxObsCoreTestMixin

TESTDIR = os.path.abspath(os.path.dirname(__file__))


class ExportRegionsTestCase(unittest.TestCase, DaxObsCoreTestMixin):
    """Tests for obscore_export_regions."""

    def setUp(self):
        self.root = makeTestTempDir(TESTDIR)

    def tearDown(self):
        removeTestTempDir(self.root)

    def _populate(self):
        butler = self.make_butler()
        self.enterContext(butler)
        butler.import_(
            filename=os.path.join(TESTDIR, "data", "hsc_gen3.yaml"),
            without_datastore=True,
        )
        return butler

    def test_visit_detector_happy_path(self):
        """Distinct visit_detector_region regions are written as STC-S."""
        self._populate()
        out = os.path.join(self.root, "regions.fits")

        obscore_export_regions(
            repo=self.root,
            dataset_type="_mock_calexp",
            collection="HSC/runs/ci_hsc",
            destination=out,
            where="",
        )

        self.assertTrue(os.path.exists(out))
        with fits.open(out) as hdul:
            self.assertEqual(len(hdul), 2)  # primary + bintable
            tbl = hdul[1]
            self.assertEqual(tbl.header["EXTNAME"], "REGIONS")
            self.assertEqual(tbl.columns.names, ["s_region"])
            self.assertEqual(tbl.header["TXTYP1"], "stc-s")
            self.assertEqual(
                tbl.header["TUTYP1"],
                "obscore:Char.SpatialAxis.Coverage.Support.Area",
            )
            self.assertEqual(tbl.header["TUCD1"], "pos.outline;obs.field")
            data = tbl.data["s_region"]
            self.assertGreater(len(data), 0)
            for s in data:
                # Every value parses as STC-S (one of the allowed leading
                # tokens). Frame keyword is always second.
                head = s.split()[0]
                self.assertIn(
                    head,
                    {"Circle", "Polygon", "Ellipse", "Union", "Intersection"},
                )
                self.assertEqual(s.split()[1], "ICRS")

    def test_patch_happy_path(self):
        """Distinct patch regions are written as STC-S."""
        self._populate()
        out = os.path.join(self.root, "regions.fits")

        obscore_export_regions(
            repo=self.root,
            dataset_type="_mock_deepCoadd",
            collection="HSC/runs/ci_hsc",
            destination=out,
            where="",
        )

        with fits.open(out) as hdul:
            data = hdul[1].data["s_region"]
            self.assertGreater(len(data), 0)
            for s in data:
                head = s.split()[0]
                self.assertIn(
                    head,
                    {"Circle", "Polygon", "Ellipse", "Union", "Intersection"},
                )

    def test_dedup_patch(self):
        """Row count equals distinct (skymap, tract, patch) keys."""
        butler = self._populate()
        out = os.path.join(self.root, "regions.fits")

        # Build the expected distinct-key count from the dataset list
        # itself by enumerating dataset refs and their dataIds.
        refs = list(butler.query_datasets("_mock_deepCoadd", collections=["HSC/runs/ci_hsc"]))
        expected_keys = {(ref.dataId["skymap"], ref.dataId["tract"], ref.dataId["patch"]) for ref in refs}
        # Sanity check: the fixture must have duplicates for this test
        # to actually exercise dedup.
        self.assertGreater(
            len(refs),
            len(expected_keys),
            "Fixture must have multiple datasets per (skymap, tract, patch) "
            "for the dedup test to be meaningful",
        )

        obscore_export_regions(
            repo=self.root,
            dataset_type="_mock_deepCoadd",
            collection="HSC/runs/ci_hsc",
            destination=out,
            where="",
        )

        with fits.open(out) as hdul:
            self.assertEqual(len(hdul[1].data), len(expected_keys))

    def test_where_filtering(self):
        """A where clause restricts the produced regions."""
        butler = self._populate()
        out_all = os.path.join(self.root, "regions_all.fits")
        out_one = os.path.join(self.root, "regions_one.fits")

        # Pick one (instrument, visit, detector) triple from the fixture
        # by enumerating the matching dataset refs.
        refs = list(butler.query_datasets("_mock_calexp", collections=["HSC/runs/ci_hsc"]))
        keys = sorted(
            {(ref.dataId["instrument"], ref.dataId["visit"], ref.dataId["detector"]) for ref in refs}
        )
        self.assertGreater(
            len(keys),
            1,
            "Need >1 distinct (visit, detector) pair in the fixture",
        )
        instrument, visit, detector = keys[0]
        where = f"instrument='{instrument}' AND visit={visit} AND detector={detector}"

        obscore_export_regions(
            repo=self.root,
            dataset_type="_mock_calexp",
            collection="HSC/runs/ci_hsc",
            destination=out_all,
            where="",
        )
        obscore_export_regions(
            repo=self.root,
            dataset_type="_mock_calexp",
            collection="HSC/runs/ci_hsc",
            destination=out_one,
            where=where,
        )

        with fits.open(out_all) as hdul_all, fits.open(out_one) as hdul_one:
            self.assertEqual(len(hdul_one[1].data), 1)
            self.assertLess(len(hdul_one[1].data), len(hdul_all[1].data))

    def test_error_no_region_dimension(self):
        """A dataset type with no region dimension is rejected."""
        butler = self._populate()
        # Register a dataset type whose dimensions have no region.
        # `instrument` alone has no region dimension.
        dt = DatasetType(
            "_mock_no_region",
            dimensions=("instrument",),
            storageClass="StructuredDataDict",
            universe=butler.dimensions,
        )
        butler.registry.registerDatasetType(dt)

        out = os.path.join(self.root, "regions.fits")
        with self.assertRaisesRegex(ValueError, "spatial region dimension"):
            obscore_export_regions(
                repo=self.root,
                dataset_type="_mock_no_region",
                collection="HSC/runs/ci_hsc",
                destination=out,
                where="",
            )

    def test_error_no_collection(self):
        """A collection is required."""
        self._populate()
        out = os.path.join(self.root, "regions.fits")

        with self.assertRaisesRegex(ValueError, "One collection"):
            obscore_export_regions(
                repo=self.root,
                dataset_type="_mock_calexp",
                collection="",
                destination=out,
                where="",
            )

    def test_empty_result(self):
        """A where clause that matches nothing produces a valid FITS file."""
        butler = self._populate()
        out = os.path.join(self.root, "regions_empty.fits")

        # Get a valid instrument from the fixture to satisfy governor
        # constraints.
        refs = list(butler.query_datasets("_mock_calexp", collections=["HSC/runs/ci_hsc"]))
        instrument = refs[0].dataId["instrument"]
        where = f"instrument='{instrument}' AND visit=-1"

        obscore_export_regions(
            repo=self.root,
            dataset_type="_mock_calexp",
            collection="HSC/runs/ci_hsc",
            destination=out,
            where=where,
        )

        with fits.open(out) as hdul:
            self.assertEqual(len(hdul[1].data), 0)
            self.assertEqual(hdul[1].columns.names, ["s_region"])
            self.assertEqual(hdul[1].header["TXTYP1"], "stc-s")

    def test_cli_smoke(self):
        """The click subcommand wiring works end-to-end."""
        self._populate()
        out = os.path.join(self.root, "regions_cli.fits")
        runner = CliRunner()
        result = runner.invoke(
            obscore_cli,
            [
                "export-regions",
                self.root,
                "_mock_calexp",
                "HSC/runs/ci_hsc",
                out,
            ],
        )
        if result.exit_code != 0:
            # Surface the exception for easier debugging.
            raise AssertionError(
                f"exit_code={result.exit_code}\noutput={result.output}\nexception={result.exception!r}"
            )
        self.assertTrue(os.path.exists(out))
        with fits.open(out) as hdul:
            self.assertGreater(len(hdul[1].data), 0)


if __name__ == "__main__":
    unittest.main()
