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

__all__ = ["DaxObsCoreTestMixin"]

from lsst.daf.butler import Butler, Config
from lsst.daf.butler.registry.obscore import DatasetTypeConfig
from lsst.dax.obscore import ExporterConfig


class DaxObsCoreTestMixin:
    """Shared code for ObsCore tests involving a butler."""

    root: str
    """Location of Butler root. Created in setUp."""

    def make_butler(self) -> Butler:
        """Return new Butler instance on each call."""
        config = Config()
        config["root"] = self.root
        config["registry", "db"] = f"sqlite:///{self.root}/gen3.sqlite3"
        butler = Butler.from_config(
            Butler.makeRepo(self.root, config=config), writeable=True, without_datastore=True
        )
        return butler

    def make_export_config(self) -> ExporterConfig:
        """Prepare configuration for exporter."""
        config = ExporterConfig(
            version=0,
            facility_name="Subaru",
            obs_collection="obs-collection",
            collections=["HSC/runs/ci_hsc"],
            use_butler_uri=False,
            dataset_types={
                "raw": DatasetTypeConfig(
                    calib_level=1,
                    dataproduct_type="image",
                    dataproduct_subtype="lsst.raw",
                    obs_id_fmt="{records[exposure].obs_id}",
                    datalink_url_fmt="http://datalink.org/{obs_id}",
                ),
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
                "r": [552.0e-9, 691.0e-9],
                "i": [691.0e-9, 818.0e-9],
                "HSC-R": [543.0e-9, 693.0e-9],
                "HSC-I": [690.0e-9, 842.0e-9],
            },
            extra_columns={"day_obs": {"template": "{records[visit].day_obs}", "type": "int"}},
        )
        return config
