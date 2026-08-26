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

import pydantic

from lsst.dax.obscore import ExporterConfig
from lsst.dax.obscore.caom_config import CaomConfig
from lsst.dax.obscore.tests import DaxObsCoreTestMixin


class CaomConfigTestCase(unittest.TestCase, DaxObsCoreTestMixin):
    """Tests of the caom configuration block."""

    def test_absent_by_default(self):
        """Configurations without a caom block are still valid."""
        config = self.make_export_config()
        self.assertIsNone(config.caom)

    def test_unknown_dataset_type_rejected(self):
        """A caom dataset type must exist in the ObsCore dataset types."""
        config = self.make_export_config()
        with self.assertRaises(pydantic.ValidationError) as cm:
            ExporterConfig.model_validate(
                config.model_dump()
                | {
                    "caom": {
                        "dataset_types": {
                            "not_a_dataset_type": {
                                "observation_id_fmt": "{visit}",
                                "product_id_fmt": "{detector}",
                            }
                        }
                    }
                }
            )
        self.assertIn("not_a_dataset_type", str(cm.exception))

    def test_algorithm_requires_derived(self):
        """An algorithm on a non-derived dataset type is a config error."""
        with self.assertRaises(pydantic.ValidationError):
            CaomConfig.model_validate(
                {
                    "dataset_types": {
                        "deep_coadd": {
                            "observation_id_fmt": "{skymap}-{tract}",
                            "product_id_fmt": "{patch}-{band}",
                            "algorithm": "some.algorithm",
                        }
                    }
                }
            )

    def test_derived_requires_algorithm(self):
        """A derived dataset type must name its algorithm."""
        with self.assertRaises(pydantic.ValidationError):
            CaomConfig.model_validate(
                {
                    "dataset_types": {
                        "deep_coadd": {
                            "observation_id_fmt": "{skymap}-{tract}",
                            "product_id_fmt": "{patch}-{band}",
                            "derived": True,
                        }
                    }
                }
            )

    def test_defaults(self):
        """Unset optional fields take the documented defaults."""
        config = CaomConfig.model_validate(
            {
                "dataset_types": {
                    "deep_coadd": {
                        "observation_id_fmt": "{skymap}-{tract}",
                        "product_id_fmt": "{patch}-{band}",
                    }
                }
            }
        )
        self.assertEqual(config.intent, "science")
        self.assertEqual(config.observation_type, "science")
        self.assertEqual(config.read_groups, [])
        self.assertIsNone(config.telescope_name)
        self.assertFalse(config.dataset_types["deep_coadd"].derived)

    def test_pixel_scale_for(self):
        """The placeholder pixel scale is read through one accessor."""
        config = CaomConfig.model_validate(
            {
                "dataset_types": {
                    "deep_coadd": {
                        "observation_id_fmt": "{skymap}-{tract}",
                        "product_id_fmt": "{patch}-{band}",
                        "s_pixel_scale": 0.2,
                    },
                    "visit_image": {
                        "observation_id_fmt": "{visit}",
                        "product_id_fmt": "{detector}",
                    },
                }
            }
        )
        self.assertEqual(config.pixel_scale_for("deep_coadd"), 0.2)
        self.assertIsNone(config.pixel_scale_for("visit_image"))
        self.assertIsNone(config.pixel_scale_for("no_such_type"))


if __name__ == "__main__":
    unittest.main()
