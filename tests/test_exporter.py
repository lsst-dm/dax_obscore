# This file is part of daf_butler.
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

import csv
import os
import unittest

import astropy.io.votable
import pyarrow
import pyarrow.parquet

from lsst.daf.butler.registry.obscore import DatasetTypeConfig
from lsst.daf.butler.registry.obscore._schema import _STATIC_COLUMNS
from lsst.daf.butler.tests.utils import makeTestTempDir, removeTestTempDir
from lsst.dax.obscore import ExporterConfig, ObscoreExporter
from lsst.dax.obscore.tests import DaxObsCoreTestMixin

TESTDIR = os.path.abspath(os.path.dirname(__file__))

# List of standard column names, some are added by a default plugin
_STANDARD_COLUMNS = tuple(col.name for col in _STATIC_COLUMNS) + ("s_dec", "s_ra", "s_fov", "s_region")


class TestCase(unittest.TestCase, DaxObsCoreTestMixin):
    """Tests of ObscoreExporter"""

    def setUp(self):
        self.root = makeTestTempDir(TESTDIR)

    def tearDown(self):
        removeTestTempDir(self.root)

    def test_schema(self):
        """Check how schema is constructed"""
        butler = self.make_butler()

        config = ExporterConfig(version=0, obs_collection="", dataset_types={}, facility_name="FACILITY")
        xprtr = ObscoreExporter(butler, config)
        self.assertCountEqual(xprtr.schema.names, _STANDARD_COLUMNS)

        # extra columns from top-level config
        config = ExporterConfig(
            version=0,
            obs_collection="",
            extra_columns={"c1": 1, "c2": "string", "c3": {"template": "{calib_level}", "type": "float"}},
            dataset_types={},
            facility_name="FACILITY",
        )
        xprtr = ObscoreExporter(butler, config)
        self.assertCountEqual(xprtr.schema.names, _STANDARD_COLUMNS + ("c1", "c2", "c3"))
        self.assertEqual(xprtr.schema.field("c1").type, pyarrow.int64())
        self.assertEqual(xprtr.schema.field("c2").type, pyarrow.string())
        self.assertEqual(xprtr.schema.field("c3").type, pyarrow.float64())

        # extra columns from per-dataset type configs
        config = ExporterConfig(
            version=0,
            obs_collection="",
            extra_columns={"c1": 1},
            dataset_types={
                "raw": DatasetTypeConfig(
                    name="raw",
                    dataproduct_type="image",
                    calib_level=1,
                    extra_columns={"c2": "string"},
                ),
                "calexp": DatasetTypeConfig(
                    dataproduct_type="image",
                    calib_level=2,
                    extra_columns={"c3": 1e10},
                ),
            },
            facility_name="FACILITY",
        )
        xprtr = ObscoreExporter(butler, config)
        self.assertCountEqual(xprtr.schema.names, _STANDARD_COLUMNS + ("c1", "c2", "c3"))
        self.assertEqual(xprtr.schema.field("c1").type, pyarrow.int64())
        self.assertEqual(xprtr.schema.field("c2").type, pyarrow.string())
        self.assertEqual(xprtr.schema.field("c3").type, pyarrow.float64())

        # Columns with the same names as in static list in configs, types
        # are not overriden.
        config = ExporterConfig(
            version=0,
            obs_collection="",
            extra_columns={"t_xel": 1e10},
            dataset_types={
                "raw": DatasetTypeConfig(
                    dataproduct_type="image",
                    calib_level=1,
                    extra_columns={"target_name": 1},
                ),
                "calexp": DatasetTypeConfig(
                    dataproduct_type="image",
                    calib_level=2,
                    extra_columns={"em_xel": "string"},
                ),
            },
            facility_name="FACILITY",
        )
        xprtr = ObscoreExporter(butler, config)
        self.assertCountEqual(xprtr.schema.names, _STANDARD_COLUMNS)
        self.assertEqual(xprtr.schema.field("t_xel").type, pyarrow.int32())
        self.assertEqual(xprtr.schema.field("target_name").type, pyarrow.string())
        self.assertEqual(xprtr.schema.field("em_xel").type, pyarrow.int32())

    def test_export_parquet(self):
        """Test Parquet export method"""
        butler = self.make_butler()
        butler.import_(filename=os.path.join(TESTDIR, "data", "hsc_gen3.yaml"), without_datastore=True)

        config = self.make_export_config()
        xprtr = ObscoreExporter(butler, config)
        output = os.path.join(self.root, "output.parquet")
        xprtr.to_parquet(output)

        table = pyarrow.parquet.read_table(output)

        def _to_python(column_name):
            """Convert table column values to Python objects."""
            for value in table.column(column_name):
                yield value.as_py()

        # Do some trivial checks
        self.assertEqual(table.num_columns, 31)
        self.assertEqual(table.num_rows, 35)
        self.assertEqual(set(_to_python("facility_name")), {"Subaru"})
        self.assertEqual(set(_to_python("obs_collection")), {"obs-collection"})
        self.assertEqual(set(_to_python("dataproduct_type")), {"image"})
        self.assertEqual(set(_to_python("dataproduct_subtype")), {"lsst.calexp", "lsst.deepCoadd"})
        self.assertEqual(set(_to_python("calib_level")), {2, 3})
        self.assertEqual(set(_to_python("instrument_name")), {"HSC", None})
        self.assertEqual(set(_to_python("em_filter_name")), {"i", "r"})
        self.assertEqual(set(_to_python("day_obs")), {20130617, 20131102, None})
        for value in _to_python("s_region"):
            self.assertTrue(value.startswith("POLYGON "))

    def test_export_csv(self):
        """Test CSV export method"""
        butler = self.make_butler()
        butler.import_(filename=os.path.join(TESTDIR, "data", "hsc_gen3.yaml"), without_datastore=True)

        # try several options for null_string
        for null_string in (None, "$NULL", ""):
            config = self.make_export_config()
            if null_string is not None:
                config.csv_null_string = null_string
                expected_null = null_string
            else:
                # default is \N
                expected_null = r"\N"

            xprtr = ObscoreExporter(butler, config)
            output = os.path.join(self.root, "output.csv")
            xprtr.to_csv(output)

            # read it back
            with open(output, newline="") as csvfile:
                reader = csv.DictReader(csvfile, delimiter=",")
                for row in reader:
                    # There should be no empty fields on output but every row
                    # is supposed to have at leas one \N value.
                    if expected_null != "":
                        self.assertFalse("" in row.values())
                    self.assertIn(expected_null, row.values())

    def test_export_votable(self):
        """Test Parquet export method"""
        butler = self.make_butler()
        butler.import_(filename=os.path.join(TESTDIR, "data", "hsc_gen3.yaml"), without_datastore=True)

        config = self.make_export_config()
        xprtr = ObscoreExporter(butler, config)
        output = os.path.join(self.root, "output.vot")
        xprtr.to_votable_file(output)

        votable = astropy.io.votable.parse(output)
        tables = list(votable.iter_tables())
        self.assertEqual(len(tables), 1)
        table0 = tables[0].array
        self.assertEqual(len(table0), 35, str(table0))

        self.assertEqual(set(table0["facility_name"]), {"Subaru"})
        self.assertEqual(set(table0["obs_collection"]), {"obs-collection"})
        self.assertEqual(set(table0["dataproduct_type"]), {"image"})
        self.assertEqual(set(table0["dataproduct_subtype"]), {"lsst.calexp", "lsst.deepCoadd"})
        self.assertEqual(set(table0["calib_level"]), {2, 3})
        self.assertEqual(set(table0["instrument_name"]), {"HSC", ""})
        self.assertEqual(set(table0["em_filter_name"]), {"i", "r"})
        self.assertEqual(set(table0["day_obs"].tolist()), {20130617, 20131102, None})
        for value in table0["s_region"]:
            self.assertTrue(value.startswith("POLYGON "))

        # Now with a limit.
        xprtr.to_votable_file(output, limit=5)
        votable = astropy.io.votable.parse(output)
        tables = list(votable.iter_tables())
        table0 = tables[0].array
        self.assertEqual(len(table0), 5, str(table0))


if __name__ == "__main__":
    unittest.main()
