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
import numbers
import sys
from collections.abc import Iterable

from lsst.daf.butler import Butler, Config

from ..config import ExporterConfig
from ..siav2 import siav2_query_from_raw

_LOG = logging.getLogger(__name__)


def obscore_siav2(
    repo: str,
    destination: str,
    config: str,
    instrument: Iterable[str],
    pos: Iterable[str],
    time: Iterable[str],
    band: Iterable[str],
    exptime: Iterable[str],
    calib: Iterable[numbers.Integral],
    maxrec: numbers.Integral | str | None,
    collections: Iterable[str],
    dataset_type: Iterable[str],
    dpsubtype: Iterable[str],
    dptype: Iterable[str],
    id: Iterable[str],
    id_list: str,
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
    instrument : `~collections.abc.Iterable` [ `str` ]
        Name of instrument to use for query.
    pos :  `~collections.abc.Iterable` [ `str` ]
        Spatial region to use for query.
    time :  `~collections.abc.Iterable` [ `str` ]
        Time or time span to use for the query, UTC MJD.
    band :  `~collections.abc.Iterable` [ `str` ]
        Wavelength range to constraint query. Units are meters.
    exptime :  `~collections.abc.Iterable` [ `str` ]
        Exposure time ranges in seconds.
    calib : `~collections.abc.Iterable` [ `int` ]
        Calibration level to select. All are selected if none specified.
    maxrec : `int` or `str` or `None`
        Maximum number of records to return. `None` is unlimited.
    collections : `~collections.abc.Iterable` [ `str` ]
        Optional collection names, if provided overrides one in ``config``.
    dataset_type : `~collections.abc.Iterable` [ `str` ]
        Names of Butler dataset types to include in query. For more general
        compatibility it is generally preferred to use dpsubtype rather than
        Butler names.
    dpsubtype : `~collections.abc.Iterable` [ `str` ]
        Data product sub types to select.
    dptype : `~collections.abc.Iterable` [ `str` ]
        The data product types to select.
    id : `~collections.abc.Iterable` [ `str` ]
        The obs_publisher_did values of explicit datasets.
    id_list : `str` or `None`
        Name of file containing multiple IDs. Will be merged with any
        explicit ID values.
    """
    config_data = Config(config)
    cfg = ExporterConfig.model_validate(config_data)

    if id_list:
        with open(id_list) as fd:
            ids = tuple(line.strip() for line in fd.readlines() if line)
        # By default `id` is a tuple. Merge with any `--id` option.
        id = tuple(id) + ids

    with Butler.from_config(repo, writeable=False) as butler:
        votable = siav2_query_from_raw(
            butler,
            cfg,
            instrument=instrument,
            pos=pos,
            time=time,
            band=band,
            exptime=exptime,
            calib=calib,
            collections=collections,
            dataset_type=dataset_type,
            dpsubtype=dpsubtype,
            dptype=dptype,
            maxrec=maxrec,
            id=id,
            query_url=" ".join(sys.argv[1:]) if sys.argv else None,
        )
    votable.to_xml(destination)
