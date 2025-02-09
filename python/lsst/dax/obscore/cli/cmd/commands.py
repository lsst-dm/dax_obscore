# This file is part of pipe_base.
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

from typing import Any

import click

from lsst.daf.butler.cli.opt import (
    collections_option,
    dataset_type_option,
    destination_argument,
    options_file_option,
    repo_argument,
    where_option,
)
from lsst.daf.butler.cli.utils import ButlerCommand, MWPath, split_commas

from ... import script


@click.group(short_help="Commands to export or update obscore data.")
def obscore() -> None:
    """Set of commands to deal with obscore data."""
    pass


@obscore.command(
    short_help="Export Butler datasets as ObsCore Data Model",
    cls=ButlerCommand,
)
@repo_argument(required=True)
@destination_argument(
    required=True,
    help="DESTINATION is the location of the output file.",
    type=MWPath(file_okay=True, dir_okay=False, writable=True),
)
@click.option(
    "--config",
    "-c",
    help="Location of the configuration file in YAML format, path or URL.",
    required=True,
)
@click.option(
    "--format",
    help="Output format, one of 'parquet', 'votable', or 'csv'; default: parquet.",
    type=click.Choice(["csv", "parquet", "votable"]),
    default="parquet",
)
@dataset_type_option(
    help=(
        "Comma-separated list of dataset types. "
        "If specified it must be a subset of dataset types defined in configuration file."
    )
)
@collections_option()
@where_option()
@options_file_option()
def export(*args: Any, **kwargs: Any) -> None:
    """Export Butler datasets as ObsCore Data Model in parquet format."""
    script.obscore_export(*args, **kwargs)


@obscore.command(
    short_help=(
        "Udpate exposure-related records (e.g. raw) in obscore table that do not have region information "
        "from matching visit-related records."
    ),
    cls=ButlerCommand,
)
@repo_argument(required=True)
@click.option("--check", help="Only print the count of records with missing regions.", is_flag=True)
@click.option("--dry-run", help="Skip actual update, but execute all other queries.", is_flag=True)
@click.option(
    "--dataproduct-type",
    help="Dataproduct type for exposure records. Default: image",
    default="image",
)
@click.option(
    "--dataproduct-subtype",
    help="Dataproduct subtype for exposure records. Default: lsst.raw",
    default="lsst.raw",
)
@click.option(
    "--instrument",
    help="If provided then updates are limited to that specific instrument.",
)
@click.option(
    "--exposure-column",
    help="Name of the column in obscore table containing exposure ID. Default: lsst_exposure",
    default="lsst_exposure",
)
@click.option(
    "--detector-column",
    help="Name of the column in obscore table containing detector ID. Default: lsst_detector",
    default="lsst_detector",
)
@click.option(
    "--region-columns",
    help=(
        "Comma-separated list of columns that all have to be NULL for updating the record. "
        "Default: s_region,s_ra,s_dec,s_fov"
    ),
    multiple=True,
    callback=split_commas,
    default=("s_region", "s_ra", "s_dec", "s_fov"),
)
def set_exposure_regions(*args: Any, **kwargs: Any) -> None:
    """Update exposure-related records (e.g. raw) in obscore table that do not
    have region information from matching visit-related records.
    """
    script.obscore_set_exposure_regions(*args, **kwargs)


@obscore.command(
    short_help=(
        "Update obscore table with the records that are missing, typically "
        "used after adding obscore support to existing repository."
    ),
    cls=ButlerCommand,
)
@repo_argument(required=True)
@click.option("--dry-run", help="Skip actual update, but execute all other queries.", is_flag=True)
def update_table(*args: Any, **kwargs: Any) -> None:
    """Update obscore table with the records that are missing, typically
    used after adding obscore support to existing repository.
    """
    script.obscore_update_table(*args, **kwargs)


@obscore.command(
    short_help="Run SIAv2 query and return results.",
    cls=ButlerCommand,
)
@repo_argument(required=True)
@destination_argument(
    required=True,
    help="DESTINATION is the location of the output file.",
    type=MWPath(file_okay=True, dir_okay=False, writable=True),
)
@click.option(
    "--config",
    "-c",
    help="Location of the configuration file in YAML format, path or URL.",
    required=True,
)
@dataset_type_option(
    help=(
        "Comma-separated list of Butler dataset types. "
        "If specified it must be a subset of dataset types defined in configuration file."
    )
)
@collections_option(
    help="Butler collections (not SIAv2 collections) to include in the query. Default is to use the "
    "collections specified in the SIAv2 configuration file."
)
@click.option(
    "--instrument",
    help="Name of instrument to use in query. If no instrument is specified all instruments are included.",
    type=str,
    multiple=True,
)
@click.option(
    "--pos",
    help="IVOA POS region to use to restrict results. CIRCLE, RANGE and POLYGON are supported.",
    type=str,
    multiple=True,
)
@click.option(
    "--time",
    help="A moment in time or a time span as a range to use to constrain the query. Uses MJD UTC.",
    type=str,
    multiple=True,
)
@click.option("--band", help="Wavelength range to constrain query. Units of meters.", type=str, multiple=True)
@click.option(
    "--exptime",
    help="Exposure time ranges in seconds.",
    type=str,
    multiple=True,
)
@click.option(
    "--calib",
    help="Calibration level of the data. Allowed values are 0, 1, 2, and 3.",
    multiple=True,
    type=int,
)
@click.option(
    "--maxrec",
    help="Maximum number of records to return. 0 means no records.",
    type=int,
)
@options_file_option()
def siav2(*args: Any, **kwargs: Any) -> None:
    """Export Butler datasets as ObsCore Data Model in VOTable format.

    For details on the SIAv2 parameters see https://www.ivoa.net/documents/SIA/

    Multiple values can be specified for a single parameter and the results are
    ORed together. All range parameters allow numbers and include -Inf and
    +Inf as options.
    """
    script.obscore_siav2(*args, **kwargs)
