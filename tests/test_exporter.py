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

import unittest

import pyarrow
from lsst.dax.obscore import DatasetTypeConfig, ExporterConfig, ObscoreExporter
from lsst.dax.obscore.obscore_exporter import _STATIC_SCHEMA


class TestCase(unittest.TestCase):
    """Tests of ObscoreExporter"""

    def test_schema(self):
        """Check how schema is constructed"""

        butler = None

        config = ExporterConfig(obs_collection="", dataset_types=[], facility_name="FACILITY")
        xprtr = ObscoreExporter(butler, config)
        self.assertEqual(xprtr.schema.names, [col[0] for col in _STATIC_SCHEMA])

        # extra columns from top-level config
        config = ExporterConfig(
            obs_collection="",
            extra_columns={"c1": 1, "c2": "string", "c3": 1e10},
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
            dataset_types=[
                DatasetTypeConfig(
                    name="raw",
                    dataproduct_type="image",
                    calib_level=1,
                    extra_columns={"c2": "string"},
                ),
                DatasetTypeConfig(
                    name="calexp",
                    dataproduct_type="image",
                    calib_level=2,
                    extra_columns={"c3": 1e10},
                ),
            ],
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
            dataset_types=[
                DatasetTypeConfig(
                    name="raw",
                    dataproduct_type="image",
                    calib_level=1,
                    extra_columns={"target_name": 1},
                ),
                DatasetTypeConfig(
                    name="calexp",
                    dataproduct_type="image",
                    calib_level=2,
                    extra_columns={"em_xel": "string"},
                ),
            ],
            facility_name="FACILITY",
        )
        xprtr = ObscoreExporter(butler, config)
        self.assertEqual(xprtr.schema.names, [col[0] for col in _STATIC_SCHEMA])
        self.assertEqual(xprtr.schema.field("t_xel").type, pyarrow.int16())
        self.assertEqual(xprtr.schema.field("target_name").type, pyarrow.string())
        self.assertEqual(xprtr.schema.field("em_xel").type, pyarrow.int16())


if __name__ == "__main__":
    unittest.main()
