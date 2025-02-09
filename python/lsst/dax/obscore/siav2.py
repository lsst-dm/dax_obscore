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

__all__ = ["SIAv2Handler", "SIAv2Parameters", "siav2_query", "siav2_query_from_raw"]

import math
import numbers
from abc import abstractmethod
from collections import defaultdict
from collections.abc import Iterable, Iterator, Sequence
from typing import Any, Self

import astropy.io.votable
import astropy.time
from pydantic import BaseModel, field_validator, model_validator

from lsst.daf.butler import Butler, DimensionGroup, Timespan
from lsst.daf.butler.pydantic_utils import SerializableRegion, SerializableTime
from lsst.sphgeom import Region, UnionRegion
from lsst.utils import inheritDoc
from lsst.utils.iteration import ensure_iterable
from lsst.utils.logging import getLogger

from .config import ExporterConfig, WhereBind
from .obscore_exporter import ObscoreExporter
from .plugins import get_siav2_handler

_LOG = getLogger(__name__)
_VALID_CALIB = frozenset({0, 1, 2, 3})


class Interval(BaseModel):
    """Representation of a simple interval."""

    start: float
    """Start of the interval."""
    end: float
    """End of the interval."""

    @model_validator(mode="after")
    def check_start_end(self) -> Self:
        """Check that start comes before end."""
        if self.start > self.end:
            raise ValueError(f"Start of interval, {self.start}, must come before end, {self.end}")
        return self

    @classmethod
    def from_string(cls, string: str) -> Self:
        """Create interval from string of form 'START END'.

        Parameters
        ----------
        string : `str`
            String representing the interval. +Inf and -Inf are allowed.
            If there is only one number the start and end interval are the
            same.

        Returns
        -------
        interval : `Interval`
            The derived interval.
        """
        interval = [float(b) for b in string.split()]
        num_values = len(interval)
        if num_values == 0:
            raise ValueError(f"No values found in supplied string {string!r}")
        if num_values > 2:
            raise ValueError(f"Found more than two values in {string!r}")
        if num_values == 1:
            interval.append(interval[0])
        return cls(start=interval[0], end=interval[1])

    def overlaps(self, other: Interval) -> bool:
        """Return whether two intervals overlap.

        Parameters
        ----------
        other : `Interval`
            Interval to check against.

        Returns
        -------
        overlaps : `bool`
            `True` if this interval overlaps the other.
        """
        return self.end >= other.start and other.end >= self.start

    def __iter__(self) -> Iterator[float]:  # type: ignore
        return iter((self.start, self.end))


class SIAv2Parameters(BaseModel):
    """Parsed versions of SIAv2 parameters."""

    pos: tuple[SerializableRegion, ...] = ()
    band: tuple[Interval, ...] = ()
    time: tuple[Timespan | SerializableTime, ...] = ()
    pol: tuple[str, ...] = ()
    fov: tuple[Interval, ...] = ()
    spatres: tuple[Interval, ...] = ()
    specrp: tuple[Interval, ...] = ()
    exptime: tuple[Interval, ...] = ()
    timeres: tuple[Interval, ...] = ()
    id: tuple[str, ...] = ()  # No validation yet.
    collection: tuple[str, ...] = ()
    facility: tuple[str, ...] = ()
    instrument: tuple[str, ...] = ()
    dptype: frozenset[str] = frozenset()
    calib: frozenset[int] = frozenset()
    target: tuple[str, ...] = ()
    maxrec: int | None = None

    @field_validator("calib")
    @classmethod
    def check_calib(cls, calib: frozenset[int]) -> frozenset[int]:
        if calib - _VALID_CALIB:
            valid = ", ".join(str(s) for s in _VALID_CALIB)
            given = ", ".join(str(s) for s in calib)
            raise ValueError(f"Calib levels can only be ({valid}) but got ({given})")
        return calib

    @field_validator("dptype")
    @classmethod
    def check_dptype(cls, dptype: frozenset[str]) -> frozenset[str]:
        valid_dptype = {"image", "cube"}
        if dptype - valid_dptype:
            raise ValueError(f"DPTYPE values can only be ({valid_dptype}) but got {dptype}")
        return dptype

    @field_validator("maxrec")
    @classmethod
    def check_maxrec(cls, maxrec: int) -> int:
        if maxrec < 0:
            raise ValueError(f"MAXREC parameter must be >= 0 but got {maxrec}")
        return maxrec

    @classmethod
    def from_siav2(
        cls,
        pos: Iterable[str] | str = (),
        band: Iterable[str] | str = (),
        time: Iterable[str] | str = (),
        pol: Iterable[str] | str = (),
        fov: Iterable[str] | str = (),
        spatres: Iterable[str] | str = (),
        specrp: Iterable[str] | str = (),
        exptime: Iterable[str] | str = (),
        timeres: Iterable[str] | str = (),
        id: Iterable[str] | str = (),
        collection: Iterable[str] | str = (),
        facility: Iterable[str] | str = (),
        instrument: Iterable[str] | str = (),
        dptype: Iterable[str] | str = (),
        calib: Iterable[numbers.Integral] | numbers.Integral = (),
        target: Iterable[str] | str = (),
        maxrec: str | numbers.Integral | None = None,
    ) -> Self:
        parsed: dict[str, Any] = {}
        if instrument:
            # Do not validate the values against the butler instruments.
            parsed["instrument"] = ensure_iterable(instrument)
        if pos:
            parsed["pos"] = [Region.from_ivoa_pos(p) for p in ensure_iterable(pos)]
        if band:
            parsed["band"] = [Interval.from_string(b) for b in ensure_iterable(band)]
        if exptime:
            parsed["exptime"] = [Interval.from_string(exp) for exp in ensure_iterable(exptime)]
        if time:
            parsed_times = []
            for tm in ensure_iterable(time):
                time_interval = Interval.from_string(tm)
                if time_interval.start == time_interval.end:
                    parsed_times.append(astropy.time.Time(time_interval.start, scale="utc", format="mjd"))
                else:
                    times: list[astropy.time.Time | None] = []
                    for t in time_interval:
                        if not math.isfinite(t):
                            # Timespan uses None to indicate unbounded.
                            times.append(None)
                        else:
                            times.append(astropy.time.Time(float(t), scale="utc", format="mjd"))
                    parsed_times.append(Timespan(times[0], times[1]))
            parsed["time"] = parsed_times
        if calib or isinstance(calib, numbers.Integral):  # 0 is allowed.
            parsed["calib"] = frozenset(int(n) for n in ensure_iterable(calib))
        if pol:
            parsed["pol"] = ensure_iterable(pol)
        if fov:
            parsed["fov"] = [Interval.from_string(val) for val in ensure_iterable(fov)]
        if spatres:
            parsed["spatres"] = [Interval.from_string(val) for val in ensure_iterable(spatres)]
        if specrp:
            parsed["specrp"] = [Interval.from_string(val) for val in ensure_iterable(specrp)]
        if timeres:
            parsed["timeres"] = [Interval.from_string(val) for val in ensure_iterable(timeres)]
        if id:
            parsed["id"] = ensure_iterable(id)
        if collection:
            parsed["collection"] = ensure_iterable(collection)
        if facility:
            parsed["facility"] = ensure_iterable(facility)
        if dptype:
            parsed["dptype"] = ensure_iterable(dptype)
        if target:
            parsed["target"] = ensure_iterable(target)
        if maxrec is not None:
            parsed["maxrec"] = int(maxrec)
        return cls.model_validate(parsed)


class SIAv2Handler:
    """Process SIAv2 query parameters and convert into Butler query
    constraints.

    Parameters
    ----------
    butler : `lsst.daf.butler.Butler`
        Butler to query.
    config : `ExporterConfig`
        ObsCore configuration.
    """

    def __init__(self, butler: Butler, config: ExporterConfig):
        self.butler = butler
        self.config = config
        self.warnings: list[str] = []

    @abstractmethod
    def get_all_instruments(self) -> list[str]:
        """Query butler for all known instruments.

        Returns
        -------
        instruments : `list` [ `str` ]
            All the instrument names known to this butler.
        """
        raise NotImplementedError()

    @abstractmethod
    def get_band_information(
        self, instruments: list[str], band_intervals: Iterable[Interval]
    ) -> dict[str, Any]:
        """Read all information from butler necessary to form a band query.

        Parameters
        ----------
        instruments : `list` [ `str` ]
            Instruments that could be involved in the band query.
        band_intervals : `~collections.abc.Iterable` of `Interval`
            The band constraints.

        Returns
        -------
        info : `dict` [ `str`, `typing.Any` ]
            Dictionary that will be provided to `from_instrument_or_band`
            when calculating the band and instrument queries.
        """
        raise NotImplementedError()

    def process_query(self, parameters: SIAv2Parameters) -> dict[str, list[WhereBind]]:
        """Process an SIAv2 query.

        Parameters
        ----------
        parameters : `SIAv2Parameters`
            Parameters to use for this query.

        Returns
        -------
        constraints : `dict` [ `str`, `list` [ `WhereBind` ]]
            The where clauses to use to find datasets, indexed by the dataset
            type name. Not all dataset types in the configuration will be
            present since some queries are incompatible with some dataset
            types.
        """
        # Reset warnings log.
        self.warnings = []
        if parameters.instrument:
            instruments = list(parameters.instrument)
        else:
            # If no explicit instrument, then all instruments could be valid.
            # Some queries like band require the instrument so ask the butler
            # for all values up front.
            instruments = self.get_all_instruments()

        band_info: dict[str, Any] = {}
        if parameters.band:
            band_info = self.get_band_information(instruments, parameters.band)

        # Loop over each dataset type calculating custom query parameters.
        dataset_type_wheres = {}
        for dataset_type_name in list(self.config.dataset_types):
            if parameters.calib:
                dataset_type_config = self.config.dataset_types[dataset_type_name]
                if dataset_type_config.calib_level not in parameters.calib:
                    continue

            wheres = []
            dataset_type = self.butler.get_dataset_type(dataset_type_name)
            dims = dataset_type.dimensions
            instrument_wheres, where = self.from_instrument_or_band(instruments, band_info, dims)
            if where is not None:
                wheres.append(where)

            if parameters.pos:
                where = self.from_pos(parameters.pos, dims)
                if not where:
                    self.warnings.append(
                        f"Can not support POS query for dataset type {dataset_type_name}. Skipping it."
                    )
                    continue
                wheres.append(where)

            if parameters.time:
                where = self.from_time(parameters.time, dims)
                if not where:
                    self.warnings.append(
                        f"Dataset type {dataset_type_name} has no timespan defined "
                        "so assuming all datasets match."
                    )
                else:
                    wheres.append(where)

            if parameters.exptime:
                exptime_wheres = self.from_exptime(parameters.exptime, dims)
                if exptime_wheres is None:
                    self.warnings.append(
                        f"Can not support EXPTIME query for dataset type {dataset_type_name}. Skipping it."
                    )
                    continue
                else:
                    wheres.append(exptime_wheres)

            # Default to an everything Where.
            where_clauses = [WhereBind(where="")]
            if instrument_wheres:
                where_clauses = [WhereBind.combine([iwhere] + wheres) for iwhere in instrument_wheres]
            elif wheres:
                where_clauses = [WhereBind.combine(wheres)]

            dataset_type_wheres[dataset_type_name] = where_clauses

        return dataset_type_wheres

    @abstractmethod
    def from_instrument_or_band(
        self, instruments: list[str], band_info: dict[str, Any], dimensions: DimensionGroup
    ) -> tuple[list[WhereBind], WhereBind | None]:
        """Convert instrument and band parameters to query information.

        Parameters
        ----------
        instruments : `list` [ `str` ]
            Instruments to include in query.
        band_info : `dict` [ `str`, `typing.Any` ]
            Matching band information. Generally obtained by calling
            `get_band_information`.
        dimensions : `lsst.daf.butler.DimensionGroup`
            The dimensions for the dataset type being queried.

        Returns
        -------
        instrument_wheres : `list` [ `WhereBind` ]
            One or more query that is instrument-specific. These will be
            treated as independent queries combined with the general query
            constraints.
        general_wheres : `WhereBind` | None
            Query constraint that is not instrument-specific.
        """
        raise NotImplementedError()

    @abstractmethod
    def from_pos(self, regions: Sequence[Region], dimensions: DimensionGroup) -> WhereBind | None:
        """Convert a region request to a butler where clause.

        Parameters
        ----------
        regions : `~collections.abc.Iterable` [ `lsst.sphgeom.Region` ]
            The region of interest.
        dimensions : `lsst.daf.butler.DimensionGroup`
            The dimensions for the dataset type being queried.

        Returns
        -------
        where : `WhereBind` or `None`
            The where clause or `None` if region is not supported by this
            dataset type.
        """
        raise NotImplementedError()

    @abstractmethod
    def from_time(
        self, ts: Iterable[astropy.time.Time | Timespan], dimensions: DimensionGroup
    ) -> WhereBind | None:
        """Convert a time span to a butler where clause.

        Parameters
        ----------
        ts : `~collections.abc.Iterable` of \
                [ `astropy.time.Time` | `lsst.daf.butler.Timespan` ]
            The times of interest.
        dimensions : `lsst.daf.butler.DimensionGroup`
            The dimensions for the dataset type being queried.

        Returns
        -------
        where : `WhereBind` or `None`
            The where clause or `None` if time is not supported by this
            dataset type.
        """
        raise NotImplementedError()

    @abstractmethod
    def from_exptime(
        self, exptime_intervals: Iterable[Interval], dimensions: DimensionGroup
    ) -> WhereBind | None:
        """Convert an exposure time interval to a butler where clause.

        Parameters
        ----------
        exptime_intervals : `~collections.abc.Iterable` [ `Interval` ]
            The exposure time interval of interest as two floating point UTC
            MJDs.
        dimensions : `lsst.daf.butler.DimensionGroup`
            The dimensions for the dataset type being queried.

        Returns
        -------
        wheres : `WhereBind` or `None`
            The where clause for the exposure times, or `None` if the dataset
            type does not understand exposure time.
        """
        raise NotImplementedError()


class SIAv2DafButlerHandler(SIAv2Handler):
    """Process SIAv2 query parameters using the daf_butler dimension universe
    and convert into Butler query constraints.
    """

    def get_all_instruments(self) -> list[str]:
        return [rec.name for rec in self.butler.query_dimension_records("instrument")]

    @inheritDoc(SIAv2Handler)
    def get_band_information(
        self, instruments: list[str], band_intervals: Iterable[Interval]
    ) -> dict[str, Any]:  # noqa: D401
        """
        Notes
        -----
        Each physical filter is checked against configuration to determine
        whether it overlaps the requested band. If it does that filter is
        stored in the returned dict with key ``filters`` where that dict maps
        instrument names to a set of matching filters. Additionally, the
        ``bands`` key is a set of bands associated with those filters
        independent of instrument.
        """  # noqa: D401
        matching_filters = defaultdict(set)
        matching_bands = set()
        with self.butler.query() as query:
            records = query.dimension_records("physical_filter")
            if instruments:
                records = records.where("instrument in (INSTRUMENTS)", bind={"INSTRUMENTS": instruments})
            for rec in records:
                if (spec_range := self.config.spectral_ranges.get(rec.name)) is not None:
                    spec_interval = Interval(start=spec_range[0], end=spec_range[1])
                    assert spec_range[0] is not None  # for mypy
                    assert spec_range[1] is not None
                    for band_interval in band_intervals:
                        if band_interval.overlaps(spec_interval):
                            matching_filters[rec.instrument].add(rec.name)
                            matching_bands.add(rec.band)
                else:
                    self.warnings.append(
                        f"Ignoring physical filter {rec.name} from instrument {rec.instrument}"
                        " since it has no defined spectral range"
                    )
        return {"filters": matching_filters, "bands": matching_bands}

    def from_instrument_or_band(
        self, instruments: list[str], band_info: dict[str, Any], dimensions: DimensionGroup
    ) -> tuple[list[WhereBind], WhereBind | None]:
        instrument_wheres = []
        general_where: WhereBind | None = None
        if "instrument" in dimensions and instruments:
            # Need separate where clauses for each instrument if
            # we also are using a physical filters constraint.
            if "physical_filter" in dimensions and "filters" in band_info:
                for i, inst in enumerate(instruments):
                    instrument_wheres.append(
                        WhereBind(
                            where=f"instrument = INST{i} AND physical_filter in (phys)",
                            bind={"phys": band_info["filters"][inst], f"INST{i}": inst},
                        )
                    )
            else:
                # Can include all instruments in single query.
                instrument_wheres.append(
                    WhereBind(where="instrument IN (INSTRUMENTS)", bind={"INSTRUMENTS": instruments})
                )
        elif "band" in dimensions and "bands" in band_info:
            # Band is not needed for an instrument query since
            # we will be using physical filters for those.
            general_where = WhereBind(where="band IN (bands)", bind={"bands": band_info["bands"]})
        return instrument_wheres, general_where

    def from_pos(self, regions: Sequence[Region], dimensions: DimensionGroup) -> WhereBind | None:
        extra_dims = set()
        region_dim = dimensions.region_dimension
        if not region_dim:
            if "exposure" in dimensions:
                # Special case for Rubin default universe.
                region_dim = "visit_detector_region"
                extra_dims = {"visit"}
            else:
                return None

        if not regions:
            # In case we are called with no regions.
            return WhereBind(where="")

        # Must OR together all regions. We can do this in sphgeom directly.
        # Butler query system itself can not OR regions.
        region = regions[0]
        if len(regions) > 1:
            for r in regions[1:]:
                region = UnionRegion(region, r)
        return WhereBind(
            where=f"{region_dim}.region OVERLAPS(region)",
            bind={"region": region},
            extra_dims=extra_dims,
        )

    def from_time(
        self, ts: Iterable[astropy.time.Time | Timespan], dimensions: DimensionGroup
    ) -> WhereBind | None:
        time_dim = dimensions.timespan_dimension
        if not time_dim:
            return None
        wheres = []
        for i, t in enumerate(ts):
            wheres.append(WhereBind(where=f"{time_dim}.timespan OVERLAPS(ts{i})", bind={f"ts{i}": t}))

        return WhereBind.combine(wheres, mode="OR")

    def from_exptime(
        self, exptime_intervals: Iterable[Interval], dimensions: DimensionGroup
    ) -> WhereBind | None:
        if "exposure" in dimensions:
            exp_dim = "exposure"
        elif "visit" in dimensions:
            exp_dim = "visit"
        else:
            return None
        wheres = []
        for exptime_interval in exptime_intervals:
            exp_wheres = []
            for cmp, value in zip((">", "<"), exptime_interval, strict=True):
                if math.isfinite(value):
                    exp_wheres.append(WhereBind(where=f"{exp_dim}.exposure_time {cmp} {value}"))
            if exp_wheres:
                wheres.append(WhereBind.combine(exp_wheres))
        if not wheres:
            # No constraint of exposure time (e.g., -Inf +Inf).
            return WhereBind()
        return WhereBind.combine(wheres, mode="OR")


def get_daf_butler_siav2_handler() -> type[SIAv2Handler]:
    """Return the SIAv2 handler specifically designed for the daf_butler
    namespace.

    Returns
    -------
    handler : `type` [ `lsst.dax.obscore.siav2.SIAv2Handler` ]
        The handler to be used for daf_butler universes.
    """
    return SIAv2DafButlerHandler


def siav2_query_from_raw(
    butler: Butler,
    config: ExporterConfig,
    *,
    pos: Iterable[str] | str = (),
    band: Iterable[str] | str = (),
    time: Iterable[str] | str = (),
    pol: Iterable[str] | str = (),
    fov: Iterable[str] | str = (),
    spatres: Iterable[str] | str = (),
    specrp: Iterable[str] | str = (),
    exptime: Iterable[str] | str = (),
    timeres: Iterable[str] | str = (),
    id: Iterable[str] | str = (),
    collection: Iterable[str] | str = (),
    facility: Iterable[str] | str = (),
    instrument: Iterable[str] | str = (),
    dptype: Iterable[str] | str = (),
    calib: Iterable[numbers.Integral] | numbers.Integral = (),
    target: Iterable[str] | str = (),
    maxrec: str | numbers.Integral | None = None,
    collections: Iterable[str] = (),
    dataset_type: Iterable[str] = (),
) -> astropy.io.votable.tree.VOTableFile:
    """Run SIAv2 query with raw parameters and return results as VOTable.

    Parameters
    ----------
    butler : `lsst.daf.butler.Butler`
        Butler repository to query.
    config : `ExporterConfig`
        Configuration for this ObsCore system.
    pos : `~collections.abc.Iterable` [ `str` ], optional
        Spatial regions to use for query. As POS region strings.
    band : `~collections.abc.Iterable` [ `str` ], optional
        Wavelength ranges to constraint query. Units are meters.
    time : `~collections.abc.Iterable` [ `str` ], optional
        Time or time spans to use for the query, UTC MJD.
    pol : `~collections.abc.Iterable` [ `str` ], optional
        Polarizations to search for.
    fov : `~collections.abc.Iterable` [ `str` ], optional
        Range of field of view to search.
    spatres : `~collections.abc.Iterable` [ `str` ], optional
        Range of spatial resolutions.
    specrp : `~collections.abc.Iterable` [ `str` ], optional
        Range of spectral resolving powers.
    exptime : `~collections.abc.Iterable` [ `str` ], optional
        Exposure time ranges in seconds.
    timeres : `~collections.abc.Iterable` [ `str` ], optional
        Range of temporal resolutions.
    id : `~collections.abc.Iterable` [ `str` ], optional
        IVO identifier of specific datasets.
    collection : `~collections.abc.Iterable` [ `str` ], optional
        ObsCore collections. It is unclear whether these are Butler collections
        or more abstract ObsCore collections that map to specific Butler
        collections.
    facility : `~collections.abc.Iterable` [ `str` ], optional
        Facility names. Currently there can only be one facility per
        ObsCore configuration of a butler so this parameter is unlikely
        to be useful.
    instrument : `~collections.abc.Iterable` [ `str` ], optional
        The names of the instrument to use to constrain the query.
    dptype : `~collections.abc.Iterable` [ `str` ], optional
        ObsCore data product type. Can be only ``image`` or ``cube``.
    calib : `~collections.abc.Iterable` [ `int` ]
        One or more calibration levels to select. Can be only 0, 1, 2, 3.
    target : `~collections.abc.Iterable` [ `str` ], optional
        Target names.
    maxrec : `str`, `numbers.Integral`, `None`, optional
        Maximum number of records to return. Unlimited if `None`. 0 means
        metadata will be returned but no records.
    collections : `~collections.abc.Iterable` [ `str` ]
        Optional collection names, if provided overrides one in ``config``.
    dataset_type : `~collections.abc.Iterable` [ `str` ]
        Names of Butler dataset types to include in query.

    Returns
    -------
    votable : `astropy.io.votable.tree.VOTableFile`
        Results of query as a VOTable.
    """
    # Parse the SIAv2 parameters
    parameters = SIAv2Parameters.from_siav2(
        instrument=instrument,
        pos=pos,
        band=band,
        time=time,
        pol=pol,
        fov=fov,
        spatres=spatres,
        specrp=specrp,
        exptime=exptime,
        timeres=timeres,
        id=id,
        collection=collection,
        facility=facility,
        calib=calib,
        dptype=dptype,
        target=target,
        maxrec=maxrec,
    )
    return siav2_query(butler, config, parameters, collections=collections, dataset_type=dataset_type)


def siav2_query(
    butler: Butler,
    config: ExporterConfig,
    parameters: SIAv2Parameters,
    *,
    collections: Iterable[str] = (),
    dataset_type: Iterable[str] = (),
) -> astropy.io.votable.tree.VOTableFile:
    """Run SIAv2 query with parsed parameters and return results as VOTable.

    Parameters
    ----------
    butler : `lsst.daf.butler.Butler`
        Butler repository to query.
    config : `ExporterConfig`
        Configuration for this ObsCore system.
    parameters : `SIAv2Parameters`
        Parsed SIAv2 parameters.
    collections : `~collections.abc.Iterable` [ `str` ]
        Optional collection names, if provided overrides one in ``config``.
    dataset_type : `~collections.abc.Iterable` [ `str` ]
        Names of dataset types to include in query.

    Returns
    -------
    votable : `astropy.io.votable.tree.VOTableFile`
        Results of query as a VOTable.
    """
    # Will modify config so copy it.
    cfg = config.model_copy(deep=True)

    if collections:
        cfg.collections = list(collections)

    if dataset_type:
        cfg.select_dataset_types(dataset_type)

    _LOG.verbose("Received parameters: %s", parameters)

    handler_type = get_siav2_handler(butler.dimensions.namespace)
    handler = handler_type(butler, cfg)
    cfg.dataset_type_constraints = handler.process_query(parameters)

    # Downselect the dataset types being queried -- a missing constraint
    # means that the dataset type is being skipped.
    cfg.select_dataset_types(cfg.dataset_type_constraints.keys())

    exporter = ObscoreExporter(butler, cfg)

    if handler.warnings:
        _LOG.warning("The query triggered the following warnings:\n%s", "\n".join(handler.warnings))

    votable = exporter.to_votable(limit=parameters.maxrec)

    if handler.warnings:
        resource = votable.resources[0]
        info = astropy.io.votable.tree.Info(name="QUERY_WARNINGS", value="OK")
        info.content = "\n".join(handler.warnings)
        resource.infos.append(info)

    return votable
