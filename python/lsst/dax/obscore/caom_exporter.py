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
import os
import re
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
# caom2.caom_util.validate_path_component rejects a component containing a
# space, slash, backslash or percent.
_INVALID_URI_COMPONENT_RE = re.compile(r"[ /\\%]")


def _sanitize_uri_component(value: str) -> str:
    """Replace characters CAOM does not allow in a URI path component.

    Parameters
    ----------
    value : `str`
        The expanded identifier.

    Returns
    -------
    sanitized : `str`
        The identifier with every space, slash, backslash and percent
        replaced by an underscore.

    Notes
    -----
    Dimension values such as skymap names legitimately contain slashes, and
    `str.format` templates offer no way to transform them, so refusing to
    export would make whole dataset types unexportable. The substitution is
    not injective, so callers must check that two different identifiers have
    not been collapsed onto one.
    """
    return _INVALID_URI_COMPONENT_RE.sub("_", value)


def _interval(lower: float, upper: float) -> Any:
    """Build a CAOM interval covering a single contiguous range.

    Parameters
    ----------
    lower : `float`
        Lower bound.
    upper : `float`
        Upper bound.

    Returns
    -------
    interval : `caom2.Interval`
        The interval, with the single sub-interval the CAOM schema
        requires.

    Notes
    -----
    The CAOM 2.4 schema makes ``samples`` mandatory on a bounds element,
    so an interval without one fails validation. Our ranges are always
    contiguous, so the sample list holds exactly one sub-interval.
    """
    return caom2.Interval(lower, upper, samples=[caom2.shape.SubInterval(lower, upper)])


def _safe_filename(observation_id: str) -> str:
    """Convert an observation identifier into a safe file name.

    Parameters
    ----------
    observation_id : `str`
        CAOM observation identifier.

    Returns
    -------
    name : `str`
        The identifier with characters that are awkward in file names
        replaced by underscores.
    """
    return re.sub(r"[^A-Za-z0-9._-]", "_", observation_id)


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

        # Records which dataset type created each plane, so that a product
        # ID collision between dataset types can be reported.
        self._plane_owners: dict[tuple[str, str], str] = {}

        # Butler URIs for the current export, keyed by dataset ID.
        self._uris: dict[Any, str] = {}

        # Source identifiers behind each sanitized identifier, so that two
        # different identifiers collapsing onto one can be reported.
        # Observation identifiers are global; product identifiers are scoped
        # to their observation.
        self._observation_id_sources: dict[str, str] = {}
        self._product_id_sources: dict[tuple[str, str], str] = {}

        # Substitutions already reported, so each is warned about once.
        self._warned_substitutions: set[tuple[str, str]] = set()

    def iter_observations(self) -> Iterator[Any]:
        """Generate CAOM Observations for the configured dataset types.

        Yields
        ------
        observation : `caom2.Observation`
            One Observation per distinct expanded ``observation_id_fmt``.
        """
        observations: dict[str, Any] = {}
        plane_data_ids: dict[tuple[str, str], Any] = {}
        state = _QueryState()
        pending: list[tuple[DatasetRef, Region | None, dict[str, Any]]] = list(
            self._obscore._iter_record_refs(state)
        )
        self._uris = self._obscore._resolve_uris([ref for ref, _, _ in pending])
        for ref, region, record in pending:
            self._add_record(observations, plane_data_ids, ref, region, record)
        self._add_auxiliary_artifacts(observations, plane_data_ids)
        yield from observations.values()

    def to_directory(self, destination: str, validate: bool = True) -> int:
        """Write one CAOM XML document per Observation.

        Parameters
        ----------
        destination : `str`
            Directory to write into. Created if it does not exist.
        validate : `bool`, optional
            If `True`, validate each document against the bundled CAOM
            schema before writing.

        Returns
        -------
        count : `int`
            Number of documents written.
        """
        os.makedirs(destination, exist_ok=True)
        writer = caom2.ObservationWriter(validate=validate)

        count = 0
        for observation in self.iter_observations():
            filename = os.path.join(destination, f"{_safe_filename(observation.observation_id)}.xml")
            with open(filename, "wb") as handle:
                writer.write(observation, handle)
            count += 1

        _LOG.info("Wrote %d CAOM observation%s to %s", count, "" if count == 1 else "s", destination)
        return count

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
        keywords.update(butler_uri=self._uris.get(ref.id, ""))
        keywords.update(record)
        return keywords

    def _add_record(
        self,
        observations: dict[str, Any],
        plane_data_ids: dict[tuple[str, str], Any],
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
        plane_data_ids : `dict` [ `tuple` [ `str`, `str` ], \
                `~lsst.daf.butler.DataCoordinate` ]
            Data ID of the primary record for each plane, keyed by
            observation and product identifier. Updated in place.
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

        observation_id = self._identifier(
            "observation_id_fmt",
            dataset_config.observation_id_fmt.format(**keywords),
            dataset_config.observation_id_fmt,
            dataset_type,
            self._observation_id_sources,
        )
        product_id = self._identifier(
            "product_id_fmt",
            dataset_config.product_id_fmt.format(**keywords),
            dataset_config.product_id_fmt,
            dataset_type,
            self._product_id_sources,
            scope=observation_id,
        )

        observation = observations.get(observation_id)
        if observation is None:
            observation = self._make_observation(observation_id, dataset_config, keywords, record)
            observations[observation_id] = observation
        else:
            self._check_observation_conflict(observation, dataset_config, keywords)

        owner_key = (observation_id, product_id)
        owner = self._plane_owners.get(owner_key)
        if owner is None:
            self._plane_owners[owner_key] = dataset_type
        elif owner != dataset_type:
            raise ValueError(
                f"Product ID {product_id!r} in observation {observation_id!r} is produced by both "
                f"{owner!r} and {dataset_type!r}. Dataset types sharing an observation must use "
                "distinct 'product_id_fmt' templates."
            )

        plane = observation.planes.get(product_id)
        if plane is None:
            plane = self._make_plane(product_id, dataset_config, ref, region, record)
            observation.planes[product_id] = plane

        plane_data_ids[owner_key] = ref.dataId

        artifact = self._make_artifact(dataset_config, keywords, product_type_name="this")
        if artifact is not None:
            plane.artifacts[artifact.uri] = artifact

    def _identifier(
        self,
        kind: str,
        raw: str,
        template: str,
        dataset_type: str,
        sources: dict[Any, str],
        scope: str | None = None,
    ) -> str:
        """Turn an expanded template into a legal CAOM identifier.

        Parameters
        ----------
        kind : `str`
            Name of the configuration key, used in messages.
        raw : `str`
            The identifier as the template produced it.
        template : `str`
            The template that produced it.
        dataset_type : `str`
            Dataset type the template belongs to.
        sources : `dict`
            Mapping from sanitized identifier to the source identifier that
            claimed it, updated in place.
        scope : `str` or `None`, optional
            Identifier this one is nested within, if any. Product
            identifiers are unique only within their observation, so they
            are keyed by ``(scope, sanitized)``.

        Returns
        -------
        sanitized : `str`
            The identifier, legal as a CAOM URI path component.

        Raises
        ------
        ValueError
            Raised if two different source identifiers sanitize to the same
            value, which would silently merge unrelated records.
        """
        sanitized = _sanitize_uri_component(raw)
        if sanitized != raw and (kind, raw) not in self._warned_substitutions:
            self._warned_substitutions.add((kind, raw))
            _LOG.warning(
                "%s %r for dataset type %s contains characters CAOM does not allow in a URI path "
                "component; exporting it as %r. Template: %s.",
                kind,
                raw,
                dataset_type,
                sanitized,
                template,
            )

        key: Any = sanitized if scope is None else (scope, sanitized)
        previous = sources.setdefault(key, raw)
        if previous != raw:
            raise ValueError(
                f"{kind} values {previous!r} and {raw!r} for dataset type {dataset_type!r} both "
                f"sanitize to {sanitized!r}, which would merge unrelated records into one. "
                f"Adjust the {kind} template {template!r} so the identifiers stay distinct."
            )
        return sanitized

    def _add_auxiliary_artifacts(
        self, observations: dict[str, Any], plane_data_ids: dict[tuple[str, str], Any]
    ) -> None:
        """Attach auxiliary dataset artifacts to the planes they belong to.

        Parameters
        ----------
        observations : `dict` [ `str`, `caom2.Observation` ]
            Observations built from the primary dataset types. Updated in
            place.
        plane_data_ids : `dict` [ `tuple` [ `str`, `str` ], \
                `~lsst.daf.butler.DataCoordinate` ]
            Data ID of the primary record for each plane, keyed by
            observation and product identifier.

        Notes
        -----
        Each auxiliary dataset type is queried once and its results indexed
        by data ID, rather than issuing a lookup for every plane.
        """
        collections = self.config.collections
        for dataset_type, config in self.caom_config.dataset_types.items():
            for auxiliary_type, product_type_name in config.auxiliary_datasets.items():
                try:
                    # with_dimension_records so auxiliary templates can use
                    # {records[...]}; limit=None so a large auxiliary
                    # dataset type is never silently truncated.
                    refs = self.butler.query_datasets(
                        auxiliary_type,
                        collections=collections,
                        with_dimension_records=True,
                        limit=None,
                        explain=False,
                    )
                except Exception as exc:
                    _LOG.warning("Could not query auxiliary dataset type %s: %s", auxiliary_type, exc)
                    continue

                by_data_id = {ref.dataId: ref for ref in refs}
                self._uris.update(self._obscore._resolve_uris(list(by_data_id.values())))

                attached = 0
                for (observation_id, product_id), data_id in plane_data_ids.items():
                    if self._plane_owners.get((observation_id, product_id)) != dataset_type:
                        continue
                    ref = by_data_id.get(data_id)
                    if ref is None:
                        continue
                    keywords = self._format_keywords(ref, {})
                    artifact = self._make_artifact(config, keywords, product_type_name)
                    if artifact is not None:
                        observations[observation_id].planes[product_id].artifacts[artifact.uri] = artifact
                        attached += 1

                _LOG.info(
                    "Attached %d auxiliary artifact%s of type %s",
                    attached,
                    "" if attached == 1 else "s",
                    auxiliary_type,
                )

    def _check_observation_conflict(
        self, observation: Any, dataset_config: CaomDatasetTypeConfig, keywords: dict[str, Any]
    ) -> None:
        """Warn if a later dataset type disagrees on observation attributes.

        Parameters
        ----------
        observation : `caom2.Observation`
            The Observation created by an earlier dataset type.
        dataset_config : `CaomDatasetTypeConfig`
            CAOM configuration for the dataset type of the current record.
        keywords : `dict` [ `str`, `~typing.Any` ]
            Template namespace for the current record.

        Notes
        -----
        The first dataset type to reach an observation ID sets the
        observation-level attributes. A later disagreement is logged and
        ignored, because there is no basis for choosing between them.
        """
        is_derived = isinstance(observation, caom2.DerivedObservation)
        if is_derived != dataset_config.derived:
            _LOG.warning(
                "Observation %s was created as %s but dataset type %s configures derived=%s; "
                "keeping the first.",
                observation.observation_id,
                type(observation).__name__,
                keywords["dataset_type"],
                dataset_config.derived,
            )

        if dataset_config.observation_type_fmt:
            observation_type = dataset_config.observation_type_fmt.format(**keywords)
            if observation_type != observation.type:
                _LOG.warning(
                    "Observation %s has type %r but dataset type %s gives %r; keeping the first.",
                    observation.observation_id,
                    observation.type,
                    keywords["dataset_type"],
                    observation_type,
                )

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
                bounds=_interval(lower, upper),
                bandpass_name=record.get("em_filter_name"),
                resolving_power=resolving_power,
                energy_bands=energy_bands,
            )

        t_min = record.get("t_min")
        t_max = record.get("t_max")
        if t_min is not None and t_max is not None:
            plane.time = caom2.Time(
                bounds=_interval(t_min, t_max),
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
