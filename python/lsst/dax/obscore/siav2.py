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

__all__ = ["SIAv2Handler", "siav2_query"]

import logging
import math
from collections import defaultdict
from collections.abc import Iterable, Iterator
from typing import Any, Self

import astropy.io.votable
import astropy.time
from lsst.daf.butler import Butler, DimensionGroup, Timespan
from lsst.daf.butler.pydantic_utils import SerializableRegion, SerializableTime
from lsst.sphgeom import Region
from pydantic import BaseModel, field_validator

from .config import ExporterConfig, WhereBind
from .obscore_exporter import ObscoreExporter

_LOG = logging.getLogger(__name__)


class Interval(BaseModel):
    """Representation of a simple interval."""

    start: float
    """Start of the interval."""
    end: float
    """End of the interval."""

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
        if len(interval) == 1:
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

    instrument: str | None = None
    pos: SerializableRegion | None = None
    time: Timespan | SerializableTime | None = None
    band: Interval | None = None
    exptime: Interval | None = None
    calib: frozenset[int] = frozenset()

    @field_validator("calib")
    @classmethod
    def check_calib(cls, calib: frozenset[int]) -> frozenset[int]:
        valid_calib = {0, 1, 2, 3}
        if calib - valid_calib:
            raise ValueError(f"Calib levels can only be ({valid_calib}) but got {calib}")
        return calib

    @classmethod
    def from_siav2(
        cls,
        instrument: str = "",
        pos: str = "",
        time: str = "",
        band: str = "",
        exptime: str = "",
        calib: Iterable[int] = (),
    ) -> Self:
        parsed: dict[str, Any] = {}
        if instrument:
            parsed["instrument"] = instrument
        if pos:
            parsed["pos"] = Region.from_ivoa_pos(pos)
        if band:
            parsed["band"] = Interval.from_string(band)
        if exptime:
            parsed["exptime"] = Interval.from_string(exptime)
        if time:
            time_interval = Interval.from_string(time)
            if time_interval.start == time_interval.end:
                parsed["time"] = astropy.time.Time(time_interval.start, scale="utc", format="mjd")
            else:
                times: list[astropy.time.Time | None] = []
                for t in time_interval:
                    if not math.isfinite(t):
                        # Timespan uses None to indicate unbounded.
                        times.append(None)
                    else:
                        times.append(astropy.time.Time(float(t), scale="utc", format="mjd"))
                parsed["time"] = Timespan(times[0], times[1])
        if calib:
            parsed["calib"] = frozenset(calib)
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

    def get_all_instruments(self) -> list[str]:
        """Query butler for all known instruments.

        Returns
        -------
        instruments : `list` [ `str` ]
            All the instrument names known to this butler.
        """
        with self.butler._query() as query:
            records = query.dimension_records("instrument")
            return [rec.name for rec in records]

    def get_band_information(self, instruments: list[str], band_interval: Interval) -> dict[str, Any]:
        """Read all information from butler necessary to form a band query.

        Parameters
        ----------
        instruments : `list` [ `str` ]
            Instruments that could be involved in the band query.
        band_interval : `Interval`
            The band constraints.

        Returns
        -------
        info : `dict` [ `str`, `typing.Any` ]
            Dictionary that will be provided to `from_instrument_or_band`
            when calculating the band and instrument queries.

        Notes
        -----
        The base class implementation assumes the default dimension universe.
        Each physical filter is checked against configuration to determine
        whether it overlaps the requested band. If it does that filter is
        stored in the returned dict with key ``filters`` where that dict maps
        instrument names to a set of matching filters. Additionally, the
        ``bands`` key is a set of bands associated with those filters
        independent of instrument.
        """
        matching_filters = defaultdict(set)
        matching_bands = set()
        with self.butler._query() as query:
            records = query.dimension_records("physical_filter")
            if instruments:
                instrs = ",".join(repr(ins) for ins in instruments)
                records = records.where(f"instrument IN ({instrs})")
            for rec in records:
                if (spec_range := self.config.spectral_ranges.get(rec.name)) is not None:
                    spec_interval = Interval(start=spec_range[0], end=spec_range[1])
                    assert spec_range[0] is not None  # for mypy
                    assert spec_range[1] is not None
                    if band_interval.overlaps(spec_interval):
                        matching_filters[rec.instrument].add(rec.name)
                        matching_bands.add(rec.band)
                else:
                    _LOG.warning(
                        "Ignoring physical filter %s since it has no defined spectral range", rec.name
                    )
        return {"filters": matching_filters, "bands": matching_bands}

    def process_query(
        self,
        instrument: str = "",
        pos: str = "",
        time: str = "",
        band: str = "",
        exptime: str = "",
        calib: Iterable[int] = (),
    ) -> dict[str, list[WhereBind]]:
        """Process an SIAv2 query.

        Parameters
        ----------
        instrument : `str`, optional
            The name of the instrument to use to constrain the query.
        pos : `str`
            Spatial region to use for query.
        time : `str`
            Time or time span to use for the query, UTC MJD.
        band : `str`
            Wavelength range to constraint query. Units are meters.
        exptime : `str`
            Exposure time ranges in seconds.
        calib : `~collections.abc.Iterable` [ `int` ]
            One or more calibration levels.

        Returns
        -------
        constraints : `dict` [ `str`, `list` [ `WhereBind` ]]
            The where clauses to use to find datasets, indexed by the dataset
            type name. Not all dataset types in the configuration will be
            present since some queries are incompatible with some dataset
            types.
        """
        # Parse the SIAv2 parameters
        parsed = SIAv2Parameters.from_siav2(
            instrument=instrument,
            pos=pos,
            band=band,
            time=time,
            exptime=exptime,
            calib=calib,
        )

        if parsed.instrument:
            instruments = [instrument]
        else:
            # If no explicit instrument, then all instruments could be valid.
            # Some queries like band require the instrument so ask the butler
            # for all values up front.
            instruments = self.get_all_instruments()

        band_info: dict[str, Any] = {}
        if parsed.band:
            band_info = self.get_band_information(instruments, parsed.band)

        # Loop over each dataset type calculating custom query parameters.
        dataset_type_wheres = {}
        for dataset_type_name in list(self.config.dataset_types):
            if parsed.calib:
                dataset_type_config = self.config.dataset_types[dataset_type_name]
                if dataset_type_config.calib_level not in parsed.calib:
                    continue

            wheres = []
            dataset_type = self.butler.get_dataset_type(dataset_type_name)
            dims = dataset_type.dimensions
            instrument_wheres, where = self.from_instrument_or_band(instruments, band_info, dims)
            if where is not None:
                wheres.append(where)

            if parsed.pos:
                where = self.from_pos(parsed.pos, dims)
                if not where:
                    _LOG.warning(
                        "Can not support POS query for dataset type %s. Skipping it.", dataset_type_name
                    )
                    continue
                wheres.append(where)

            if parsed.time:
                where = self.from_time(parsed.time, dims)
                if not where:
                    _LOG.warning(
                        "Dataset type %s has no timespan defined so assuming all datasets match.",
                        dataset_type_name,
                    )
                else:
                    wheres.append(where)

            if parsed.exptime:
                exptime_wheres = self.from_exptime(parsed.exptime, dims)
                if exptime_wheres is None:
                    _LOG.warning(
                        "Can not support EXPTIME query for dataset type %s. Skipping it.", dataset_type_name
                    )
                    continue
                else:
                    wheres.extend(exptime_wheres)

            if instrument_wheres:
                where_clauses = [WhereBind.combine([iwhere] + wheres) for iwhere in instrument_wheres]
            else:
                where_clauses = [WhereBind.combine(wheres)]

            dataset_type_wheres[dataset_type_name] = where_clauses

        return dataset_type_wheres

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
        instrument_wheres = []
        general_where: WhereBind | None = None
        if "instrument" in dimensions and instruments:
            # Need separate where clauses for each instrument if
            # we also are using a physical filters constraint.
            if "physical_filter" in dimensions and "filters" in band_info:
                for inst in instruments:
                    instrument_wheres.append(
                        WhereBind(
                            where=f"instrument = {inst!r} AND physical_filter in (phys)",
                            bind={"phys": band_info["filters"][inst]},
                        )
                    )
            else:
                # Can include all instruments in query, although
                # binding does not work.
                instrs = ",".join(repr(inst) for inst in instruments)
                instrument_wheres.append(WhereBind(where=f"instrument IN ({instrs})"))
        elif "band" in dimensions and "bands" in band_info:
            # Band is not needed for an instrument query since
            # we will be using physical filters for those.
            general_where = WhereBind(where="band IN (bands)", bind={"bands": band_info["bands"]})
        return instrument_wheres, general_where

    def from_pos(self, region: Region, dimensions: DimensionGroup) -> WhereBind | None:
        """Convert a region request to a butler where clause.

        Parameters
        ----------
        region : `lsst.sphgeom.Region`
            The region of interest.
        dimensions : `lsst.daf.butler.DimensionGroup`
            The dimensions for the dataset type being queried.

        Returns
        -------
        where : `WhereBind` or `None`
            The where clause or `None` if region is not supported by this
            dataset type.
        """
        extra_dims = set()
        region_dim_element = dimensions.region_dim
        if not region_dim_element:
            if "exposure" in dimensions:
                # Special case for Rubin default universe.
                region_dim = "visit_detector_region"
                extra_dims = {"visit"}
            else:
                return None
        else:
            region_dim = region_dim_element.name
        return WhereBind(
            where=f"{region_dim}.region OVERLAPS(region)",
            bind={"region": region},
            extra_dims=extra_dims,
        )

    def from_time(self, ts: astropy.time.Time | Timespan, dimensions: DimensionGroup) -> WhereBind | None:
        """Convert a time span to a butler where clause.

        Parameters
        ----------
        ts : `astropy.time.Time` or `lsst.daf.butler.Timespan`
            The time of interest.
        dimensions : `lsst.daf.butler.DimensionGroup`
            The dimensions for the dataset type being queried.

        Returns
        -------
        where : `WhereBind` or `None`
            The where clause or `None` if time is not supported by this
            dataset type.
        """
        time_dim = dimensions.timespan_dim
        if not time_dim:
            return None
        return WhereBind(where=f"{time_dim}.timespan OVERLAPS(ts)", bind={"ts": ts})

    def from_exptime(self, exptime_interval: Interval, dimensions: DimensionGroup) -> list[WhereBind] | None:
        """Convert an exposure time interval to a butler where clause.

        Parameters
        ----------
        exptime_interval : `Interval`
            The exposure time interval of interest as two floating point UTC
            MJDs.
        dimensions : `lsst.daf.butler.DimensionGroup`
            The dimensions for the dataset type being queried.

        Returns
        -------
        wheres : `list` [ `WhereBind` ] or `None`
            The where clauses for the exposure times, or `None` if the dataset
            type does not understand exposure time.
        """
        if "exposure" in dimensions:
            exp_dim = "exposure"
        elif "visit" in dimensions:
            exp_dim = "visit"
        else:
            return None
        wheres = []
        for cmp, value in zip((">", "<"), exptime_interval, strict=True):
            if math.isfinite(value):
                wheres.append(WhereBind(where=f"{exp_dim}.exposure_time {cmp} {value}"))
        return wheres


def siav2_query(
    butler: Butler,
    config: ExporterConfig,
    instrument: str = "",
    pos: str = "",
    time: str = "",
    band: str = "",
    exptime: str = "",
    calib: Iterable[int] = (),
    collections: Iterable[str] = (),
    dataset_type: Iterable[str] = (),
) -> astropy.io.votable.tree.VOTableFile:
    """Export Butler datasets as ObsCore Data Model in parquet format.

    Parameters
    ----------
    butler : `lsst.daf.butler.Butler`
        Butler repository to query.
    config : `ExporterConfig`
        Configuration for this ObsCore system.
    instrument : `str`
        Name of instrument to use for query.
    pos : `str`
        Spatial region to use for query.
    time : `str`
        Time or time span to use for the query, UTC MJD.
    band : `str`
        Wavelength range to constraint query. Units are meters.
    exptime : `str`
        Exposure time ranges in seconds.
    calib : `~collections.abc.Iterable` [ `str` ]
       Calibration levels to select.
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

    handler = SIAv2Handler(butler, cfg)
    cfg.dataset_type_constraints = handler.process_query(instrument, pos, time, band, exptime, calib)

    # Downselect the dataset types being queried -- a missing constraint
    # means that the dataset type is being skipped.
    cfg.select_dataset_types(cfg.dataset_type_constraints.keys())

    exporter = ObscoreExporter(butler, cfg)
    return exporter.to_votable()
