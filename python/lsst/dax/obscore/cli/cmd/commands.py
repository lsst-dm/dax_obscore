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
from lsst.daf.butler.cli.utils import ButlerCommand, MWPath

from ... import script


@click.command(
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
    help="Output format, one of 'parquet' or 'csv'; default: parquet.",
    type=click.Choice(["csv", "parquet"]),
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
def obscore_export(*args: Any, **kwargs: Any) -> None:
    """Export Butler datasets as ObsCore Data Model in parquet format."""
    script.obscore_export(*args, **kwargs)
