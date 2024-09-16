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

__all__ = ["obscore_update_table"]

import logging
import re
import time
from collections.abc import Iterator
from typing import Any

from lsst.daf.butler import Butler, CollectionType, DatasetType
from lsst.daf.butler.registry.interfaces import CollectionRecord
from lsst.daf.butler.registry.obscore import (
    ConfigCollectionType,
    ObsCoreLiveTableManager,
    ObsCoreManagerConfig,
)
from lsst.daf.butler.registry.sql_registry import SqlRegistry
from lsst.utils import iteration

_LOG = logging.getLogger(__name__)


def obscore_update_table(
    repo: str,
    dry_run: bool,
) -> None:
    """Update obscore table with the records that are missing, typically
    used after adding obscore support to existing repository.

    Parameters
    ----------
    repo : `str`
        URI to the butler repository.
    dry_run : `bool`
        If `True` then print the records that will be added, but do not
        actually add them.
    """
    butler = Butler.from_config(repo, writeable=True)

    # There is no client API for updating obscore table, so we need to access
    # internals of the Registry and obscore manager.
    # Have to use non-public Registry interface.
    registry = butler._registry  # type: ignore
    assert isinstance(registry, SqlRegistry), "Registry must be SqlRegistry"
    manager = registry._managers.obscore
    assert manager is not None, "Registry is not configured for obscore support"
    assert isinstance(manager, ObsCoreLiveTableManager), f"Unexpected type of obscore manager {type(manager)}"
    config = manager.config

    for collection_record, dataset_type in _collections(registry, config):
        start_time = time.time()
        refs = registry.queryDatasets(dataset_type, collections=collection_record.name).expanded()
        if dry_run:
            for ref in refs:
                _LOG.info("Will be adding dataset %s", ref)
        else:
            count = 0
            if collection_record.type is CollectionType.RUN:
                # Limit record number in single insert.
                for refs_chunk in iteration.chunk_iterable(refs):
                    count += manager.add_datasets(refs_chunk)
            elif collection_record.type is CollectionType.TAGGED:
                for refs_chunk in iteration.chunk_iterable(refs):
                    count += manager.associate(refs_chunk, collection_record)
            else:
                raise ValueError(f"Unexpected collection type: {collection_record.type}")
            end_time = time.time()
            _LOG.info(
                "Added %s records for dataset type %r and collection %r in %.0f seconds",
                count,
                dataset_type.name,
                collection_record.name,
                end_time - start_time,
            )


def _collections(
    registry: SqlRegistry, config: ObsCoreManagerConfig
) -> Iterator[tuple[CollectionRecord, DatasetType]]:
    """Generate list of collections and dataset types to search for dataset
    references.

    Yields
    ------
    collection_record : `CollectionRecord`
        Record for a collection to search.
    dataset_type : `DatasetType`
        Dataset type to search.
    """
    dataset_types = list(registry.queryDatasetTypes(config.dataset_types.keys()))
    _LOG.info("Found these dataset types in registry: %s", [ds.name for ds in dataset_types])

    collections: Any
    if config.collection_type is ConfigCollectionType.RUN:
        collection_type = CollectionType.RUN
        if config.collections is None:
            collections = ...
        else:
            collections = [re.compile(collection) for collection in config.collections]
    elif config.collection_type is ConfigCollectionType.TAGGED:
        collection_type = CollectionType.TAGGED
        collections = config.collections
    else:
        raise ValueError(f"Unexpected collection type: {config.collection_type}")

    for dataset_type in dataset_types:
        for collection in registry.queryCollections(
            collections, datasetType=dataset_type, collectionTypes=collection_type
        ):
            collection_record = registry._managers.collections.find(collection)
            yield collection_record, dataset_type
