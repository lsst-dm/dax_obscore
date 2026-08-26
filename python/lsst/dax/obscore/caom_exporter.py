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

__all__ = ["CaomExporter"]

import datetime
from collections.abc import Iterator
from typing import Any

import astropy.units as u
from astropy.coordinates import EarthLocation

from lsst.daf.butler import Butler, DatasetRef
from lsst.sphgeom import Region
from lsst.utils.logging import getLogger

from ._caom_shapes import arcsec_to_degrees, caom2, metres, region_to_caom_shape, require_caom2
from .caom_config import CaomDatasetTypeConfig
from .config import ExporterConfig
from .obscore_exporter import ObscoreExporter, _QueryState

_LOG = getLogger(__name__)

# CAOM identifiers become components of a "caom:<collection>/<id>" URI, and
# caom2.caom_util.validate_path_component rejects a component containing any
# of these characters.
_INVALID_URI_COMPONENT_CHARS = (" ", "/", "\\", "%")


def _validate_uri_component(kind: str, value: str, template: str, dataset_type: str) -> str:
    """Check an expanded identifier is legal in a CAOM URI.

    Parameters
    ----------
    kind : `str`
        Name of the configuration key, used in the error message.
    value : `str`
        The expanded identifier.
    template : `str`
        The template that produced it.
    dataset_type : `str`
        Dataset type the template belongs to.

    Returns
    -------
    value : `str`
        The identifier, unchanged.

    Raises
    ------
    ValueError
        Raised if the identifier contains a character CAOM does not allow.

    Notes
    -----
    The identifier is not rewritten to make it legal, because it is the
    identifier CADC ingests; the configuration must produce a legal one.
    """
    bad = [char for char in _INVALID_URI_COMPONENT_CHARS if char in value]
    if bad:
        raise ValueError(
            f"{kind} {value!r} for dataset type {dataset_type!r} contains {bad}, which CAOM does not "
            f"allow in a URI path component (space, slash, backslash and percent are forbidden). "
            f"Adjust the {kind} template {template!r}."
        )
    return value


class CaomExporter:
    """Export Butler datasets as CAOM Observations.

    Parameters
    ----------
    butler : `lsst.daf.butler.Butler`
        Data butler.
    config : `lsst.dax.obscore.ExporterConfig`
        Exporter configuration. Must have a ``caom`` block.

    Raises
    ------
    ValueError
        Raised if the configuration has no ``caom`` block.
    ImportError
        Raised if the optional ``caom2`` dependency is not installed.
    """

    def __init__(self, butler: Butler, config: ExporterConfig):
        require_caom2()
        if config.caom is None:
            raise ValueError("Configuration has no 'caom' section; cannot export CAOM.")

        self.butler = butler
        self.config = config
        self.caom_config = config.caom

        # Only export dataset types with CAOM configuration.
        config.select_dataset_types(self.caom_config.dataset_types)
        self._obscore = ObscoreExporter(butler, config)

    def iter_observations(self) -> Iterator[Any]:
        """Generate CAOM Observations for the configured dataset types.

        Yields
        ------
        observation : `caom2.Observation`
            One Observation per distinct expanded ``observation_id_fmt``.
        """
        observations: dict[str, Any] = {}
        state = _QueryState()
        for ref, region, record in self._obscore._iter_record_refs(state):
            self._add_record(observations, ref, region, record)
        yield from observations.values()

    def _format_keywords(self, ref: DatasetRef, record: dict[str, Any]) -> dict[str, Any]:
        """Build the namespace used to expand configuration templates.

        Parameters
        ----------
        ref : `~lsst.daf.butler.DatasetRef`
            Reference to the dataset.
        record : `dict` [ `str`, `~typing.Any` ]
            The completed ObsCore record.

        Returns
        -------
        keywords : `dict` [ `str`, `~typing.Any` ]
            Values available to every template, comprising the data ID
            mapping, its dimension records, the dataset identity, and every
            ObsCore column.
        """
        keywords: dict[str, Any] = {"records": ref.dataId.records}
        keywords.update(ref.dataId.mapping)
        keywords.update(id=ref.id, run=ref.run, dataset_type=ref.datasetType.name)
        keywords.update(record)
        return keywords

    def _add_record(
        self,
        observations: dict[str, Any],
        ref: DatasetRef,
        region: Region | None,
        record: dict[str, Any],
    ) -> None:
        """Fold one ObsCore record into the accumulating Observations.

        Parameters
        ----------
        observations : `dict` [ `str`, `caom2.Observation` ]
            Observations accumulated so far, keyed by observation ID.
            Updated in place.
        ref : `~lsst.daf.butler.DatasetRef`
            Reference to the dataset.
        region : `~lsst.sphgeom.Region` or `None`
            Spatial region for the dataset.
        record : `dict` [ `str`, `~typing.Any` ]
            The completed ObsCore record.
        """
        dataset_type = ref.datasetType.name
        dataset_config = self.caom_config.dataset_types[dataset_type]
        keywords = self._format_keywords(ref, record)

        observation_id = _validate_uri_component(
            "observation_id_fmt",
            dataset_config.observation_id_fmt.format(**keywords),
            dataset_config.observation_id_fmt,
            dataset_type,
        )
        product_id = _validate_uri_component(
            "product_id_fmt",
            dataset_config.product_id_fmt.format(**keywords),
            dataset_config.product_id_fmt,
            dataset_type,
        )

        observation = observations.get(observation_id)
        if observation is None:
            observation = self._make_observation(observation_id, dataset_config, keywords, record)
            observations[observation_id] = observation

        plane = observation.planes.get(product_id)
        if plane is None:
            plane = self._make_plane(product_id, dataset_config, ref, region, record)
            observation.planes[product_id] = plane

        artifact = self._make_artifact(dataset_config, keywords, product_type_name="this")
        if artifact is not None:
            plane.artifacts[artifact.uri] = artifact

    def _make_observation(
        self,
        observation_id: str,
        dataset_config: CaomDatasetTypeConfig,
        keywords: dict[str, Any],
        record: dict[str, Any],
    ) -> Any:
        """Create a CAOM Observation.

        Parameters
        ----------
        observation_id : `str`
            Expanded observation identifier.
        dataset_config : `CaomDatasetTypeConfig`
            CAOM configuration for the dataset type of the first record
            encountered for this observation.
        keywords : `dict` [ `str`, `~typing.Any` ]
            Template namespace for this record.
        record : `dict` [ `str`, `~typing.Any` ]
            The completed ObsCore record.

        Returns
        -------
        observation : `caom2.Observation`
            A new Observation with no planes.
        """
        config = self.caom_config
        origin = self.config.origin

        observation_type = (
            dataset_config.observation_type_fmt.format(**keywords)
            if dataset_config.observation_type_fmt
            else config.observation_type
        )
        intent_name = (
            dataset_config.intent_fmt.format(**keywords) if dataset_config.intent_fmt else config.intent
        )
        intent = caom2.ObservationIntentType(intent_name.lower())

        target = None
        if record.get("target_name"):
            target = caom2.Target(name=record["target_name"])

        instrument = None
        if record.get("instrument_name"):
            instrument = caom2.Instrument(record["instrument_name"])

        proposal = None
        if config.proposal_id is not None:
            proposal = caom2.Proposal(
                id=config.proposal_id,
                project=self.config.obs_collection,
                title=origin.title if origin is not None else None,
            )

        meta_release = None
        if origin is not None:
            meta_release = datetime.datetime.combine(
                origin.publication_date, datetime.time(), tzinfo=datetime.UTC
            )

        kwargs: dict[str, Any] = {
            "collection": self.config.obs_collection,
            "observation_id": observation_id,
            "intent": intent,
            "type": observation_type,
            "proposal": proposal,
            "telescope": self._make_telescope(record),
            "instrument": instrument,
            "target": target,
            "meta_release": meta_release,
        }
        if dataset_config.derived:
            return caom2.DerivedObservation(algorithm=caom2.Algorithm(dataset_config.algorithm), **kwargs)
        return caom2.SimpleObservation(**kwargs)

    def _make_telescope(self, record: dict[str, Any]) -> Any:
        """Create a CAOM Telescope with a geocentric position.

        Parameters
        ----------
        record : `dict` [ `str`, `~typing.Any` ]
            The completed ObsCore record, used for the facility name.

        Returns
        -------
        telescope : `caom2.Telescope`
            The telescope, with geocentric coordinates when they can be
            resolved.
        """
        facility = record.get("facility_name") or self.config.facility_name
        name = self.caom_config.telescope_name or facility

        location = self.caom_config.geo_location
        if location is not None:
            x, y, z = location
        else:
            try:
                site = EarthLocation.of_site(facility)
            except Exception:
                _LOG.warning(
                    "Could not resolve a geocentric position for facility %r; "
                    "set caom.geo_location to supply one.",
                    facility,
                )
                return caom2.Telescope(name=name)
            x, y, z = (metres(value) for value in site.geocentric)

        return caom2.Telescope(name=name, geo_location_x=x, geo_location_y=y, geo_location_z=z)

    def _make_plane(
        self,
        product_id: str,
        dataset_config: CaomDatasetTypeConfig,
        ref: DatasetRef,
        region: Region | None,
        record: dict[str, Any],
    ) -> Any:
        """Create a CAOM Plane with its position, energy and time.

        Parameters
        ----------
        product_id : `str`
            Expanded product identifier.
        dataset_config : `CaomDatasetTypeConfig`
            CAOM configuration for the dataset type.
        ref : `~lsst.daf.butler.DatasetRef`
            Reference to the dataset.
        region : `~lsst.sphgeom.Region` or `None`
            Spatial region for the dataset.
        record : `dict` [ `str`, `~typing.Any` ]
            The completed ObsCore record.

        Returns
        -------
        plane : `caom2.Plane`
            A new plane with no artifacts.
        """
        config = self.caom_config
        origin = self.config.origin

        provenance = None
        if config.provenance is not None:
            provenance = caom2.Provenance(
                name=config.provenance.name,
                version=config.provenance.version,
                project=config.provenance.project,
                producer=origin.publisher if origin is not None else None,
                run_id=ref.run,
                reference=origin.reference_url if origin is not None else None,
            )

        observable = None
        if record.get("o_ucd"):
            observable = caom2.Observable(record["o_ucd"])

        read_groups = caom2.caom_util.URISet()
        for uri in config.read_groups:
            read_groups.add(uri)

        plane = caom2.Plane(
            product_id=product_id,
            data_product_type=caom2.DataProductType(record["dataproduct_type"]),
            calibration_level=caom2.CalibrationLevel(record["calib_level"]),
            provenance=provenance,
            observable=observable,
            data_read_groups=read_groups,
        )

        bounds = region_to_caom_shape(region)
        dimension = None
        if record.get("s_xel1") is not None and record.get("s_xel2") is not None:
            dimension = caom2.Dimension2D(record["s_xel1"], record["s_xel2"])
        sample_size = arcsec_to_degrees(config.pixel_scale_for(ref.datasetType.name))
        if bounds is not None or dimension is not None or sample_size is not None:
            plane.position = caom2.Position(
                bounds=bounds,
                dimension=dimension,
                resolution=record.get("s_resolution"),
                sample_size=sample_size,
            )

        em_min = record.get("em_min")
        em_max = record.get("em_max")
        if em_min is not None and em_max is not None:
            # ObsCore em_min and em_max are already metres, as CAOM
            # expects, but the conversion is stated rather than assumed.
            lower = metres(em_min * u.m)
            upper = metres(em_max * u.m)
            resolving_power = ((upper + lower) / 2.0) / (upper - lower) if upper > lower else None
            # em_band is deprecated in CAOM 2.4 and removed in 2.5, so the
            # band is supplied through energy_bands instead.
            energy_bands = None
            if config.em_band:
                energy_bands = caom2.caom_util.TypedSet(caom2.EnergyBand, caom2.EnergyBand[config.em_band])
            plane.energy = caom2.Energy(
                bounds=caom2.Interval(lower, upper),
                bandpass_name=record.get("em_filter_name"),
                resolving_power=resolving_power,
                energy_bands=energy_bands,
            )

        t_min = record.get("t_min")
        t_max = record.get("t_max")
        if t_min is not None and t_max is not None:
            plane.time = caom2.Time(
                bounds=caom2.Interval(t_min, t_max),
                exposure=record.get("t_exptime"),
            )

        return plane

    def _make_artifact(
        self, dataset_config: CaomDatasetTypeConfig, keywords: dict[str, Any], product_type_name: str
    ) -> Any | None:
        """Create a CAOM Artifact for one dataset.

        Parameters
        ----------
        dataset_config : `CaomDatasetTypeConfig`
            CAOM configuration for the dataset type.
        keywords : `dict` [ `str`, `~typing.Any` ]
            Template namespace for this record.
        product_type_name : `str`
            CAOM product type name. Must be a member value of
            `caom2.ProductType`, which has no ``science`` term: the primary
            artifact of a plane uses ``this`` and extras use ``auxiliary``.

        Returns
        -------
        artifact : `caom2.Artifact` or `None`
            The artifact, or `None` if no URI could be determined.
        """
        uri_fmt = self.caom_config.artifact_uri_fmt
        if uri_fmt is None:
            _LOG.warning(
                "No caom.artifact_uri_fmt configured and no Butler URI available for %s; skipping artifact.",
                keywords["id"],
            )
            return None
        uri = uri_fmt.format(**keywords)

        return caom2.Artifact(
            uri=uri,
            product_type=caom2.ProductType(product_type_name),
            release_type=caom2.ReleaseType.DATA,
            content_type=dataset_config.content_type or self.caom_config.content_type,
        )
