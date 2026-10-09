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

__all__ = ["obscore_export_caom"]

from collections.abc import Iterable

from lsst.daf.butler import Butler, Config

from ..caom_exporter import CaomExporter
from ..config import ExporterConfig, WhereBind


def obscore_export_caom(
    repo: str,
    destination: str,
    config: str,
    where: str | None,
    collections: Iterable[str],
    dataset_type: Iterable[str],
) -> None:
    """Export Butler datasets as CAOM observations.

    Parameters
    ----------
    repo : `str`
        URI to the butler repository.
    destination : `str`
        Directory to write the CAOM XML documents into.
    config : `str`
        Location of the configuration file.
    where : `str` or `None`
        Optional user expression, if provided overrides one in ``config``.
    collections : `~collections.abc.Iterable` [ `str` ]
        Optional collection names, if provided overrides those in
        ``config``.
    dataset_type : `~collections.abc.Iterable` [ `str` ]
        Names of dataset types to export. Must be a subset of the dataset
        types configured in the ``caom`` section.
    """
    config_data = Config(config)
    cfg = ExporterConfig.model_validate(config_data)
    if cfg.caom is None:
        raise ValueError(f"Configuration {config} has no 'caom' section; cannot export CAOM.")
    if where:
        cfg.where = WhereBind(where=where)
    if collections:
        cfg.collections = list(collections)
    if dataset_type:
        requested = set(dataset_type)
        unknown = requested - set(cfg.caom.dataset_types)
        if unknown:
            raise ValueError(f"Dataset types {sorted(unknown)} have no 'caom' configuration.")
        cfg.caom.dataset_types = {
            name: value for name, value in cfg.caom.dataset_types.items() if name in requested
        }

    with Butler.from_config(repo, writeable=False) as butler:
        CaomExporter(butler, cfg).to_directory(destination)
