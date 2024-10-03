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

__all__ = ["obscore_export"]

from collections.abc import Iterable

from lsst.daf.butler import Butler, Config

from ..config import ExporterConfig, WhereBind
from ..obscore_exporter import ObscoreExporter


def obscore_export(
    repo: str,
    destination: str,
    config: str,
    format: str,
    where: str | None,
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
    where : `str`
        Optional user expression, if provided overrides one in ``config``.
    collections : `~collections.abc.Iterable` [ `str` ]
        Optional collection names, if provided overrides one in ``config``.
    dataset_type : `~collections.abc.Iterable` [ `str` ]
        Names of dataset types to export.
    """
    butler = Butler.from_config(repo, writeable=False)

    config_data = Config(config)
    cfg = ExporterConfig.model_validate(config_data)
    if where:
        cfg.where = WhereBind(where=where)
    if collections:
        cfg.collections = list(collections)
    if dataset_type:
        cfg.select_dataset_types(dataset_type)

    exporter = ObscoreExporter(butler, cfg)
    if format == "parquet":
        exporter.to_parquet(destination)
    elif format == "csv":
        exporter.to_csv(destination)
    elif format == "votable":
        exporter.to_votable_file(destination)
    else:
        raise ValueError(f"Unexpected output format {format:r}")
