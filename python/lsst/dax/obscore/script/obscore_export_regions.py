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

__all__ = ["obscore_export_regions"]

import logging

import astropy.io.fits as fits
import numpy as np

from lsst.daf.butler import Butler

_LOG = logging.getLogger(__name__)


def obscore_export_regions(
    repo: str,
    dataset_type: str,
    collection: str,
    destination: str,
    where: str,
) -> None:
    """Export distinct spatial regions for the datasets matching a query
    as a FITS binary table of IVOA STC-S strings.

    Parameters
    ----------
    repo : `str`
        URI to the butler repository.
    dataset_type : `str`
        Name of the butler dataset type whose regions should be exported.
    collection : `str`
        Collection to search.
    destination : `str`
        Path to the output FITS file.
    where : `str`
        Butler query ``where`` expression. Empty string means no constraint.
    """
    if not collection:
        raise ValueError("One collection must be specified.")

    with Butler.from_config(repo, writeable=False) as butler:
        dt = butler.get_dataset_type(dataset_type)
        region_dim = dt.dimensions.region_dimension
        if region_dim is None:
            raise ValueError(
                f"Dataset type {dataset_type!r} has no single unambiguous spatial region dimension."
            )
        region_key = f"{region_dim}.region"
        region_group = butler.dimensions[region_dim].minimal_group

        with butler.query() as query:
            q = query.join_dataset_search(dt.name, collections=[collection])
            if where:
                q = q.where(where)
            result = q.general(region_group, region_key)
            stcs_strings = [row[region_key].to_ivoa_stcs() for row in result]

    _write_fits(stcs_strings, destination, dataset_type, collection, where)
    _LOG.info("Wrote %d region(s) to %s", len(stcs_strings), destination)


def _write_fits(
    stcs_strings: list[str],
    destination: str,
    dataset_type: str,
    collection: str,
    where: str,
) -> None:
    """Write the STC-S strings as a single-extension FITS binary table."""
    if stcs_strings:
        max_len = max(len(s) for s in stcs_strings)
    else:
        max_len = 1
    arr = np.array(stcs_strings, dtype=f"U{max_len}")
    col = fits.Column(name="s_region", format=f"{max_len}A", array=arr)
    hdu = fits.BinTableHDU.from_columns([col])
    hdu.header["EXTNAME"] = "REGIONS"
    hdu.header["TUTYP1"] = "obscore:Char.SpatialAxis.Coverage.Support.Area"
    hdu.header["TUCD1"] = "pos.outline;obs.field"
    hdu.header["TXTYP1"] = "stc-s"
    hdu.header["HISTORY"] = (
        f'export-regions dataset_type={dataset_type} collection={collection} where="{where}"'
    )
    hdul = fits.HDUList([fits.PrimaryHDU(), hdu])
    hdul.writeto(destination, overwrite=True)
