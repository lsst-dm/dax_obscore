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

__all__ = ["process_siav2_parameters", "siav2_query"]

import logging
import math
from collections import defaultdict
from collections.abc import Iterable

import astropy.io.votable
import astropy.time
from lsst.daf.butler import Butler, Timespan
from lsst.sphgeom import Region

from .config import ExporterConfig, WhereBind
from .obscore_exporter import ObscoreExporter

_LOG = logging.getLogger(__name__)


def _overlaps(start1: float, end1: float, start2: float, end2: float) -> bool:
    """Return whether range (start1, end1) overlaps with (start2, end2)

    Parameters
    ----------
    start1 : `float`
        Start of first range.
    end1 : `float`
        End of first range.
    start2 : `float`
        Start of second range.
    end2 : `float`
        End of second range.

    Returns
    -------
    overlaps : `bool`
        `True` if range 1 overlaps range 2.
    """
    return end1 >= start2 and end2 >= start1


def process_siav2_parameters(
    butler: Butler,
    cfg: ExporterConfig,
    instrument: str,
    pos: str,
    time: str,
    band: str,
    exptime: str,
) -> dict[str, list[WhereBind]]:
    """Process SIAv2 parameters and calculate user where expressions
    for each dataset type.

    Parameters
    ----------
    butler : `lsst.daf.butler.Butler`
        The butler to query for additional information.
    cfg : `ExporterConfig`
        Export config associated with this query.
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

    Returns
    -------
    wheres : `dict` [ `str`, `list` [ `WhereBind` ]]
        A dictionary with keys of the dataset type names and one or more
        where clauses to apply when querying.
    """
    # Some queries require instruments and if no instrument is specified
    # all instruments must be queried.
    if instrument:
        instruments = [instrument]
    else:
        with butler._query() as query:
            records = query.dimension_records("instrument")
            instruments = [rec.name for rec in records]

    siav2_bind = {}  # Global binds available to queries.
    if pos:
        siav2_bind["region"] = Region.from_ivoa_pos(pos)
    if time:
        components = [float(t) for t in time.split()]
        if len(components) == 1:
            siav2_bind["ts"] = astropy.time.Time(components[0], scale="utc", format="mjd")
        elif len(components) == 2:
            times: list[astropy.time.Time | None] = []
            for t in components:
                if not math.isfinite(t):
                    # Timespan uses None to indicate unbounded.
                    times.append(None)
                else:
                    times.append(astropy.time.Time(float(t), scale="utc", format="mjd"))
            siav2_bind["ts"] = Timespan(times[0], times[1])
        else:
            raise ValueError("Too many times in TIME field.")
    if band:
        # BAND is a wavelength in m.
        ranges = [float(b) for b in band.split()]
        if len(ranges) == 1:
            ranges.append(ranges[0])
        # Configuration maps from physical_filter to wavelength.
        # This works for instrument/visit/exposure dimensions.
        # For coadds there is only band (ugrizy) and no instrument.
        matching_filters = defaultdict(set)
        matching_bands = set()
        with butler._query() as query:
            records = query.dimension_records("physical_filter")
            if instrument:
                records = records.where("instrument = instr", bind={"instr": instrument})
            for rec in records:
                if rec.name in cfg.spectral_ranges:
                    spec_range = cfg.spectral_ranges[rec.name]
                    assert spec_range[0] is not None  # for mypy
                    assert spec_range[1] is not None
                    if _overlaps(ranges[0], ranges[1], spec_range[0], spec_range[1]):
                        matching_filters[rec.instrument].add(rec.name)
                        matching_bands.add(rec.band)
                else:
                    _LOG.warning(
                        "Ignoring physical filter %s since it has no defined spectral range", rec.name
                    )
        siav2_bind["filters"] = matching_filters
        siav2_bind["bands"] = matching_bands

    dataset_type_wheres = {}
    for dataset_type_name in list(cfg.dataset_types):
        wheres = []
        dataset_type = butler.get_dataset_type(dataset_type_name)
        dims = dataset_type.dimensions
        instrument_wheres = []
        if "instrument" in dims and instruments:
            # Need separate where clauses for each instrument if
            # we also are using a physical filters constraint.
            if "physical_filter" in dims and "filters" in siav2_bind:
                for inst in instruments:
                    instrument_wheres.append(
                        WhereBind(
                            where=f"instrument = {inst!r} AND physical_filter in (phys)",
                            bind={"phys": siav2_bind["filters"][inst]},
                        )
                    )
            else:
                # Can include all instruments in query, although
                # binding does not work.
                instrs = ",".join(repr(inst) for inst in instruments)
                instrument_wheres.append(WhereBind(where=f"instrument IN ({instrs})", bind={}))
        elif "band" in dims and "bands" in siav2_bind:
            # Band is not needed for an instrument query since
            # we will be using physical filters for those.
            wheres.append(WhereBind(where="band IN (bands)", bind={"bands": siav2_bind["bands"]}))
        if pos:
            extra_dims = set()
            region_dim_element = dims.region_dim
            if not region_dim_element:
                if "exposure" in dims:
                    # Special case for Rubin default universe.
                    region_dim = "visit_detector_region"
                    extra_dims = {"visit"}
                else:
                    _LOG.warning("Can not support POS query for dataset type %s", dataset_type_name)
                    del cfg.dataset_types[dataset_type_name]
                    continue
            else:
                region_dim = region_dim_element.name
            wheres.append(
                WhereBind(
                    where=f"{region_dim}.region OVERLAPS(region)",
                    bind={"region": siav2_bind["region"]},
                    extra_dims=extra_dims,
                )
            )
        if time:
            time_dim = dims.timespan_dim
            if not time_dim:
                _LOG.warning(
                    "Dataset type %s has no timespan defined so assuming all datasets match.",
                    dataset_type_name
                )
            else:
                wheres.append(
                    WhereBind(where=f"{time_dim}.timespan OVERLAPS(ts)", bind={"ts": siav2_bind["ts"]})
                )
        if exptime:
            if "exposure" in dims:
                exp_dim = "exposure"
            elif "visit" in dims:
                exp_dim = "visit"
            else:
                _LOG.warning("Can not support EXPTIME query for dataset type %s", dataset_type_name)
                del cfg.dataset_types[dataset_type_name]
                continue
            ranges = [float(e) for e in exptime.split()]
            if math.isfinite(ranges[0]):
                wheres.append(WhereBind(where=f"{exp_dim}.exposure_time > {ranges[0]}", bind={}))
            if math.isfinite(ranges[1]):
                wheres.append(WhereBind(where=f"{exp_dim}.exposure_time < {ranges[1]}", bind={}))

        if instrument_wheres:
            where_clauses = [WhereBind.combine([iwhere] + wheres) for iwhere in instrument_wheres]
        else:
            where_clauses = [WhereBind.combine(wheres)]

        dataset_type_wheres[dataset_type_name] = where_clauses

    return dataset_type_wheres


def siav2_query(
    butler: Butler,
    config: ExporterConfig,
    instrument: str = "",
    pos: str = "",
    time: str = "",
    band: str = "",
    exptime: str = "",
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

    cfg.dataset_type_constraints = process_siav2_parameters(butler, cfg, instrument, pos, time, band, exptime)

    exporter = ObscoreExporter(butler, cfg)
    return exporter.to_votable()
