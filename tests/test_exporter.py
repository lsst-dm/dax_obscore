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

import os
import unittest

import pyarrow
import pyarrow.parquet
from lsst.daf.butler import Butler, Config
from lsst.daf.butler.tests import DatastoreMock
from lsst.daf.butler.tests.utils import makeTestTempDir, removeTestTempDir
from lsst.dax.obscore import DatasetTypeConfig, ExporterConfig, ObscoreExporter
from lsst.dax.obscore.obscore_exporter import _STATIC_SCHEMA

TESTDIR = os.path.abspath(os.path.dirname(__file__))


class TestCase(unittest.TestCase):
    """Tests of ObscoreExporter"""

    def setUp(self):
        self.root = makeTestTempDir(TESTDIR)

    def tearDown(self):
        removeTestTempDir(self.root)

    def make_butler(self) -> Butler:
        """Return new Butler instance on each call."""
        config = Config()
        config["root"] = self.root
        config["registry", "db"] = f"sqlite:///{self.root}/gen3.sqlite3"
        butler = Butler(Butler.makeRepo(self.root, config=config), writeable=True)
        DatastoreMock.apply(butler)
        return butler

    def test_schema(self):
        """Check how schema is constructed"""

        butler = None

        config = ExporterConfig(obs_collection="", dataset_types=[], facility_name="FACILITY")
        xprtr = ObscoreExporter(butler, config)
        self.assertEqual(xprtr.schema.names, [col[0] for col in _STATIC_SCHEMA])

        # extra columns from top-level config
        config = ExporterConfig(
            obs_collection="",
            extra_columns={"c1": 1, "c2": "string", "c3": {"template": "{calib_level}", "type": "float"}},
            dataset_types=[],
            facility_name="FACILITY",
        )
        xprtr = ObscoreExporter(butler, config)
        self.assertEqual(
            xprtr.schema.names,
            [col[0] for col in _STATIC_SCHEMA] + ["c1", "c2", "c3"],
        )
        self.assertEqual(xprtr.schema.field("c1").type, pyarrow.int64())
        self.assertEqual(xprtr.schema.field("c2").type, pyarrow.string())
        self.assertEqual(xprtr.schema.field("c3").type, pyarrow.float64())

        # extra columns from per-dataset type configs
        config = ExporterConfig(
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
        self.assertEqual(
            xprtr.schema.names,
            [col[0] for col in _STATIC_SCHEMA] + ["c1", "c2", "c3"],
        )
        self.assertEqual(xprtr.schema.field("c1").type, pyarrow.int64())
        self.assertEqual(xprtr.schema.field("c2").type, pyarrow.string())
        self.assertEqual(xprtr.schema.field("c3").type, pyarrow.float64())

        # Columns with the same names as in static list in configs, types
        # are not overriden.
        config = ExporterConfig(
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
        self.assertEqual(xprtr.schema.names, [col[0] for col in _STATIC_SCHEMA])
        self.assertEqual(xprtr.schema.field("t_xel").type, pyarrow.int16())
        self.assertEqual(xprtr.schema.field("target_name").type, pyarrow.string())
        self.assertEqual(xprtr.schema.field("em_xel").type, pyarrow.int16())

    def test_export(self):
        """Test export method"""
        butler = self.make_butler()
        butler.import_(filename=os.path.join(TESTDIR, "data", "hsc_gen3.yaml"))

        config = ExporterConfig(
            facility_name="Subaru",
            obs_collection="obs-collection",
            collections=["HSC/runs/ci_hsc"],
            use_butler_uri=False,
            dataset_types={
                "_mock_calexp": DatasetTypeConfig(
                    calib_level=2,
                    dataproduct_type="image",
                    dataproduct_subtype="lsst.calexp",
                    obs_id_fmt="{records[visit].name}",
                    datalink_url_fmt="http://datalink.org/{obs_id}",
                ),
                "_mock_deepCoadd": DatasetTypeConfig(
                    calib_level=3,
                    dataproduct_type="image",
                    dataproduct_subtype="lsst.deepCoadd",
                    obs_id_fmt="{skymap}-{tract}-{patch}",
                    datalink_url_fmt="http://datalink.org/{id}",
                ),
            },
            spectral_ranges={
                "r": [552.0, 691.0],
                "i": [691.0, 818.0],
            },
            extra_columns={"day_obs": {"template": "{records[visit].day_obs}", "type": "int"}},
        )
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


if __name__ == "__main__":
    unittest.main()
