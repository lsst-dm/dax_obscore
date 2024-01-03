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

from .. import ExporterConfig, ObscoreExporter


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
        Output format, 'csv' or 'parquet'
    where : `str`
        Optional user expression, if provided overrides one in ``config``.
    collections : `iterable` [ `str` ]
        Optional collection names, if provided overrides one in ``config``.
    """
    butler = Butler(repo, writeable=False)

    config_data = Config(config)
    cfg = ExporterConfig.parse_obj(config_data)
    if where:
        cfg.where = where
    if collections:
        cfg.collections = list(collections)
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
    else:
        raise ValueError(f"Unexpected output format {format:r}")
