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

__all__ = ["CaomConfig", "CaomDatasetTypeConfig", "CaomProvenanceConfig"]

from pydantic import BaseModel, Field, model_validator


class CaomProvenanceConfig(BaseModel):
    """Provenance values that apply to every plane in the export."""

    name: str
    """Collection-specific common name of the process that produced the
    data, for example ``LSST Science Pipelines``. This names the software,
    not the dataset."""

    version: str | None = None
    """Version of that software, for example ``v29.0``. This is a fixed
    release-level value; the per-execution value is the Butler RUN
    collection, which is reported as ``Provenance.runID``."""

    project: str | None = None
    """Project the process belongs to, for example ``DRP``."""


class CaomDatasetTypeConfig(BaseModel):
    """CAOM configuration for a single dataset type."""

    observation_id_fmt: str
    """Template for the CAOM ``Observation.observationID``. Dataset types
    whose templates expand to the same value contribute planes to the same
    Observation."""

    product_id_fmt: str
    """Template for the CAOM ``Plane.productID``, unique within an
    Observation."""

    algorithm: str | None = None
    """Algorithm name for a derived Observation. Only valid when
    ``derived`` is true."""

    derived: bool = False
    """If true, emit a ``DerivedObservation`` rather than a
    ``SimpleObservation``."""

    content_type: str | None = None
    """MIME type of the files, overriding the top-level default."""

    s_pixel_scale: float | None = None
    """Pixel scale in arcsec, used for ``Position.sampleSize``.

    This is a placeholder. The value belongs beside ``s_xel`` in the
    ObsCore ``DatasetTypeConfig``, and moves there once daf_butler
    supports it. Read it through `CaomConfig.pixel_scale_for` so that the
    migration touches one function.
    """

    observation_type_fmt: str | None = None
    """Template for ``Observation.type``, overriding the top-level
    literal. Use this for exposure-based dataset types, where the value
    is carried on the exposure dimension record."""

    intent_fmt: str | None = None
    """Template for ``Observation.intent``, overriding the top-level
    literal. Must expand to ``science`` or ``calibration``."""

    auxiliary_datasets: dict[str, str] = Field(default_factory=dict)
    """Dataset types contributing auxiliary artifacts to each plane.

    Keys are dataset type names. Values must be member values of
    `caom2.ProductType`, which does not include ``science``; use
    ``auxiliary``, ``weight``, ``preview`` and so on. The plane's primary
    artifact always uses ``this``.
    """

    @model_validator(mode="after")
    def _check_algorithm(self) -> CaomDatasetTypeConfig:
        """Reject an algorithm on a non-derived dataset type.

        Returns
        -------
        config : `CaomDatasetTypeConfig`
            This configuration, unchanged.
        """
        if self.algorithm is not None and not self.derived:
            raise ValueError(
                "'algorithm' is only meaningful when 'derived' is true; "
                "caom2 forces a SimpleObservation to use the algorithm name 'exposure'."
            )
        if self.derived and self.algorithm is None:
            raise ValueError("'derived' dataset types must specify an 'algorithm'.")
        return self


class CaomConfig(BaseModel):
    """CAOM extensions to an ObsCore exporter configuration.

    Notes
    -----
    Values already present in the ObsCore configuration are not repeated
    here. ``obs_collection`` supplies the CAOM collection,
    ``facility_name`` the telescope name, and the ``origin`` block the
    proposal title, producer, provenance reference and metadata release
    date.
    """

    telescope_name: str | None = None
    """CAOM ``Telescope.name``. Defaults to the ObsCore facility name."""

    geo_location: tuple[float, float, float] | None = None
    """Geocentric telescope position in metres, overriding the lookup by
    facility name. Needed only for facilities astropy does not know."""

    proposal_id: str | None = None
    """CAOM ``Proposal.id``."""

    em_band: str | None = None
    """CAOM ``EnergyBand`` name, for example ``OPTICAL``."""

    provenance: CaomProvenanceConfig | None = None
    """Provenance values shared by every plane."""

    read_groups: list[str] = Field(default_factory=list)
    """Group URIs controlling access to the data."""

    intent: str = "science"
    """Default ``Observation.intent``."""

    observation_type: str = "science"
    """Default ``Observation.type``.

    This matches the vocabulary of the ``exposure`` dimension's
    ``observation_type`` field, so that literal and templated values agree.
    """

    content_type: str | None = None
    """Default MIME type of the exported files."""

    artifact_uri_fmt: str | None = None
    """Template for ``Artifact.uri``. ``{butler_uri}`` is available when a
    datastore is present. If unset, the resolved Butler URI is used."""

    dataset_types: dict[str, CaomDatasetTypeConfig]
    """Per-dataset-type configuration. Only these dataset types are
    exported."""

    def pixel_scale_for(self, dataset_type: str) -> float | None:
        """Return the pixel scale in arcsec for a dataset type.

        Parameters
        ----------
        dataset_type : `str`
            Name of the dataset type.

        Returns
        -------
        pixel_scale : `float` or `None`
            Pixel scale in arcsec, or `None` if not configured.

        Notes
        -----
        This is the single point of access for the placeholder
        ``s_pixel_scale`` key. When the value moves to the ObsCore
        ``DatasetTypeConfig``, only this method changes.
        """
        config = self.dataset_types.get(dataset_type)
        return config.s_pixel_scale if config is not None else None
