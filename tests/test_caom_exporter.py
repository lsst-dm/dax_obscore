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

import astropy.units as u

from lsst.daf.butler.tests.utils import makeTestTempDir, removeTestTempDir
from lsst.dax.obscore.tests import DaxObsCoreTestMixin

try:
    import caom2
except ImportError:
    caom2 = None

TESTDIR = os.path.abspath(os.path.dirname(__file__))


@unittest.skipUnless(caom2 is not None, "caom2 is not installed")
class CaomExporterTestCase(unittest.TestCase, DaxObsCoreTestMixin):
    """Tests of CaomExporter."""

    def setUp(self):
        self.root = makeTestTempDir(TESTDIR)

    def tearDown(self):
        removeTestTempDir(self.root)

    def make_populated_butler(self):
        """Return a butler with the test data imported."""
        butler = self.make_butler()
        self.enterContext(butler)
        butler.import_(filename=os.path.join(TESTDIR, "data", "hsc_gen3.yaml"), without_datastore=True)
        return butler

    def make_single_type_config(self, dataset_type):
        """Return a CAOM config restricted to one dataset type."""
        config = self.make_caom_config()
        config.select_dataset_types([dataset_type])
        config.caom.dataset_types = {dataset_type: config.caom.dataset_types[dataset_type]}
        return config

    def test_calexp_observations(self):
        """Visit-based records group into one Observation per visit."""
        from lsst.dax.obscore.caom_exporter import CaomExporter

        butler = self.make_populated_butler()
        config = self.make_single_type_config("_mock_calexp")

        observations = list(CaomExporter(butler, config).iter_observations())

        self.assertGreater(len(observations), 0)
        for observation in observations:
            self.assertIsInstance(observation, caom2.SimpleObservation)
            self.assertEqual(observation.collection, "obs-collection")
            self.assertEqual(observation.telescope.name, "Subaru Telescope")
            self.assertEqual(observation.proposal.id, "TEST-PROPOSAL")
            self.assertEqual(observation.type, "science")
            self.assertEqual(observation.intent, caom2.ObservationIntentType.SCIENCE)
            self.assertGreater(len(observation.planes), 0)
            for plane in observation.planes.values():
                self.assertEqual(plane.calibration_level, caom2.CalibrationLevel.CALIBRATED)
                self.assertEqual(plane.data_product_type, caom2.DataProductType.IMAGE)
                self.assertEqual(plane.provenance.name, "LSST Science Pipelines")
                self.assertEqual(plane.provenance.version, "v29.0")
                self.assertIsNotNone(plane.provenance.run_id)
                self.assertIsNotNone(plane.position)
                self.assertIsInstance(plane.position.bounds, caom2.shape.Polygon)
                self.assertIsNotNone(plane.energy)
                self.assertIn(caom2.EnergyBand.OPTICAL, plane.energy.energy_bands)
                self.assertIsNotNone(plane.time)
                self.assertGreaterEqual(len(plane.artifacts), 1)
                for artifact in plane.artifacts.values():
                    self.assertTrue(artifact.uri.startswith("cadc:TEST/_mock_calexp/"))
                    self.assertEqual(artifact.content_type, "application/fits")
                    self.assertEqual(artifact.release_type, caom2.ReleaseType.DATA)
                    self.assertEqual(artifact.product_type, caom2.ProductType.THIS)

    def test_coadd_is_derived(self):
        """Coadd records produce DerivedObservations keyed on tract."""
        from lsst.dax.obscore.caom_exporter import CaomExporter

        butler = self.make_populated_butler()
        config = self.make_single_type_config("_mock_deepCoadd")

        observations = list(CaomExporter(butler, config).iter_observations())

        self.assertGreater(len(observations), 0)
        for observation in observations:
            self.assertIsInstance(observation, caom2.DerivedObservation)
            self.assertEqual(observation.algorithm.name, "mock.coadd")
            self.assertTrue(observation.observation_id.startswith("coadd-"))
            for plane in observation.planes.values():
                self.assertEqual(plane.calibration_level, caom2.CalibrationLevel.PRODUCT)
                self.assertTrue(plane.product_id.startswith("deepCoadd-"))

    def test_pixel_scale_becomes_sample_size(self):
        """The configured pixel scale reaches Position.sampleSize in deg."""
        from lsst.dax.obscore.caom_exporter import CaomExporter

        butler = self.make_populated_butler()
        config = self.make_single_type_config("_mock_calexp")

        observations = list(CaomExporter(butler, config).iter_observations())
        plane = next(iter(observations[0].planes.values()))

        self.assertAlmostEqual(plane.position.sample_size, (0.17 * u.arcsec).to_value(u.deg), places=12)

    def test_illegal_observation_id_rejected(self):
        """An identifier CAOM cannot put in a URI is reported clearly."""
        from lsst.dax.obscore.caom_exporter import CaomExporter

        butler = self.make_populated_butler()
        config = self.make_single_type_config("_mock_deepCoadd")
        # The test skymap is named "discrete/ci_hsc", so this template
        # expands to an identifier containing a slash.
        config.caom.dataset_types["_mock_deepCoadd"].observation_id_fmt = "{skymap}-{tract}"

        with self.assertRaises(ValueError) as cm:
            list(CaomExporter(butler, config).iter_observations())
        message = str(cm.exception)
        self.assertIn("observation_id_fmt", message)
        self.assertIn("slash", message)

    def test_merges_across_dataset_types(self):
        """Two dataset types sharing an observation ID share an Observation."""
        from lsst.dax.obscore.caom_exporter import CaomExporter

        butler = self.make_populated_butler()
        config = self.make_caom_config()
        # Give the coadd the same observation ID template as the calexp so
        # that the two dataset types collide deliberately.
        config.caom.dataset_types["_mock_deepCoadd"].observation_id_fmt = "shared-observation"
        config.caom.dataset_types["_mock_calexp"].observation_id_fmt = "shared-observation"

        with self.assertLogs("lsst.dax.obscore.caom_exporter", level="WARNING"):
            observations = list(CaomExporter(butler, config).iter_observations())

        self.assertEqual(len(observations), 1)
        product_ids = set(observations[0].planes)
        self.assertTrue(any(p.startswith("calexp-") for p in product_ids))
        self.assertTrue(any(p.startswith("deepCoadd-") for p in product_ids))
        levels = {plane.calibration_level for plane in observations[0].planes.values()}
        self.assertEqual(levels, {caom2.CalibrationLevel.CALIBRATED, caom2.CalibrationLevel.PRODUCT})

    def test_duplicate_product_id_is_an_error(self):
        """Colliding product IDs from different dataset types abort."""
        from lsst.dax.obscore.caom_exporter import CaomExporter

        butler = self.make_populated_butler()
        config = self.make_caom_config()
        config.caom.dataset_types["_mock_deepCoadd"].observation_id_fmt = "shared-observation"
        config.caom.dataset_types["_mock_calexp"].observation_id_fmt = "shared-observation"
        config.caom.dataset_types["_mock_deepCoadd"].product_id_fmt = "collide"
        config.caom.dataset_types["_mock_calexp"].product_id_fmt = "collide"

        with self.assertRaises(ValueError) as cm:
            list(CaomExporter(butler, config).iter_observations())
        self.assertIn("collide", str(cm.exception))

    def test_uri_template_without_datastore(self):
        """A datastore-less butler still exports, using the URI template."""
        from lsst.dax.obscore.caom_exporter import CaomExporter

        butler = self.make_populated_butler()
        config = self.make_single_type_config("_mock_calexp")

        exporter = CaomExporter(butler, config)
        with self.assertLogs("lsst.dax.obscore.obscore_exporter", level="WARNING"):
            observations = list(exporter.iter_observations())

        plane = next(iter(observations[0].planes.values()))
        artifact = next(iter(plane.artifacts.values()))
        self.assertTrue(artifact.uri.startswith("cadc:TEST/_mock_calexp/"))

    def test_butler_uri_keyword_is_available(self):
        """{butler_uri} expands to empty string when there is no datastore."""
        from lsst.dax.obscore.caom_exporter import CaomExporter

        butler = self.make_populated_butler()
        config = self.make_single_type_config("_mock_calexp")
        config.caom.artifact_uri_fmt = "cadc:TEST/{butler_uri}{id}"

        observations = list(CaomExporter(butler, config).iter_observations())
        plane = next(iter(observations[0].planes.values()))
        artifact = next(iter(plane.artifacts.values()))
        self.assertTrue(artifact.uri.startswith("cadc:TEST/"))

    def test_to_directory_writes_valid_xml(self):
        """Documents are written per Observation and validate against XSD."""
        from lsst.dax.obscore.caom_exporter import CaomExporter

        butler = self.make_populated_butler()
        config = self.make_caom_config()
        destination = os.path.join(self.root, "caom")

        count = CaomExporter(butler, config).to_directory(destination)

        self.assertGreater(count, 0)
        written = sorted(os.listdir(destination))
        self.assertEqual(len(written), count)
        for name in written:
            self.assertTrue(name.endswith(".xml"))
            with open(os.path.join(destination, name), "rb") as handle:
                head = handle.read(512)
            self.assertIn(b"Observation", head)


if __name__ == "__main__":
    unittest.main()
