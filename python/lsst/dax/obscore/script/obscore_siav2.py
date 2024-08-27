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

__all__ = ["obscore_siav2"]

import logging
from collections.abc import Iterable

from lsst.daf.butler import Butler, Config

from .. import ExporterConfig, ObscoreExporter

_LOG = logging.getLogger(__name__)


def _overlaps(start1: float, end1: float, start2: float, end2: float) -> bool:
    """Return whether range (start1, end1) overlaps with (start2, end2)?"""
    return end1 >= start2 and end2 >= start1


def obscore_siav2(
    repo: str,
    destination: str,
    config: str,
    format: str,
    instrument: str,
    pos: str,
    time: str,
    band: str,
    collections: Iterable[str],
    dataset_type: Iterable[str],
) -> None:
    """Export Butler datasets as ObsCore Data Model in parquet format.

    Parameters
    ----------
    repo : `str`
        URI to the butler repository.
    destination : `str`
        Location of the output file.
    config : `str`
        Location of the configuration file.
    format : `str`
        Output format, 'csv' or 'parquet'.
    instrument : `str`
        Name of instrument to use for query.
    pos : `str`
        Spatial region to use for query.
    time : `str`
        Time or time span to use for the query, UTC MJD.
    band : `str`
        Wavelength range to constraint query. Units are meters.
    collections : `~collections.abc.Iterable` [ `str` ]
        Optional collection names, if provided overrides one in ``config``.
    dataset_type : `~collections.abc.Iterable` [ `str` ]
        Names of dataset types to include in query.
    """
    butler = Butler.from_config(repo, writeable=False)

    config_data = Config(config)
    cfg = ExporterConfig.model_validate(config_data)
    if pos:
        cfg.siav2["POS"] = pos
    if time:
        cfg.siav2["TIME"] = time
    if instrument:
        cfg.siav2["INSTRUMENT"] = instrument
    if collections:
        cfg.collections = list(collections)

    if band:
        # BAND is a wavelength in m.
        ranges = [float(b) for b in band.split()]
        if len(ranges) == 1:
            ranges.append(ranges[0])
        # Configuration maps from physical_filter to wavelength.
        # This works for instrument/visit/exposure dimensions.
        # For coadds there is only band (ugrizy) and no instrument.
        matching_bands = {}
        with butler._query() as query:
            records = query.dimension_records("physical_filter")
            for rec in records:
                if instrument and rec.instrument != instrument:
                    continue
                if rec.name in cfg.spectral_ranges:
                    if _overlaps(*ranges, *cfg.spectral_ranges[rec.name]):
                        # Ignore instrument for now.
                        matching_bands[rec.name] = rec.band
                else:
                    _LOG.warning(
                        "Ignoring physical filter %s since it has no defined spectral range", rec.name
                    )
        cfg.siav2["filters"] = list(matching_bands.keys())
        cfg.siav2["bands"] = set(matching_bands.values())

    if dataset_type:
        dataset_type_set = set(dataset_type)
        # Check that configuration has all requested dataset types.
        if not dataset_type_set.issubset(cfg.dataset_types):
            extras = dataset_type_set - set(cfg.dataset_types)
            raise ValueError(f"Dataset types {extras} are not defined in configuration file.")
        # Remove dataset types that are not needed.
        cfg.dataset_types = {
            key: value for key, value in cfg.dataset_types.items() if key in dataset_type_set
        }

    exporter = ObscoreExporter(butler, cfg)
    if format == "parquet":
        exporter.to_parquet(destination)
    elif format == "csv":
        exporter.to_csv(destination)
    elif format == "votable":
        exporter.to_votable(destination)
    else:
        raise ValueError(f"Unexpected output format {format:r}")
