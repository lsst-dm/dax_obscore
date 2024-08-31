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

from ..config import ExporterConfig
from ..siav2 import siav2_query

_LOG = logging.getLogger(__name__)


def obscore_siav2(
    repo: str,
    destination: str,
    config: str,
    format: str,
    instrument: str,
    pos: str,
    time: str,
    band: str,
    exptime: str,
    calib: Iterable[int],
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
    exptime : `str`
        Exposure time ranges in seconds.
    calib : `~collections.abc.Iterable` [ `int` ]
        Calibration level to select. All are selected if none specified.
    collections : `~collections.abc.Iterable` [ `str` ]
        Optional collection names, if provided overrides one in ``config``.
    dataset_type : `~collections.abc.Iterable` [ `str` ]
        Names of dataset types to include in query.
    """
    butler = Butler.from_config(repo, writeable=False)

    config_data = Config(config)
    cfg = ExporterConfig.model_validate(config_data)

    votable = siav2_query(butler, cfg, instrument, pos, time, band, exptime, calib, collections, dataset_type)
    votable.to_xml(destination)
