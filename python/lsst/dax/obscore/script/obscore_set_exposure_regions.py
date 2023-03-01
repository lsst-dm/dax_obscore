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

__all__ = ["obscore_set_exposure_regions"]

import logging
from collections import defaultdict
from collections.abc import Collection, MutableMapping
from typing import Any

import sqlalchemy
from lsst.daf.butler import Butler, DatasetId
from lsst.daf.butler.registries.sql import SqlRegistry
from lsst.daf.butler.registry.obscore import ObsCoreLiveTableManager, Record
from lsst.sphgeom import Region
from lsst.utils.iteration import chunk_iterable

_LOG = logging.getLogger(__name__)


def obscore_set_exposure_regions(
    repo: str,
    check: bool,
    dry_run: bool,
    dataproduct_type: str,
    dataproduct_subtype: str,
    instrument: str | None,
    exposure_column: str,
    detector_column: str,
    region_columns: Collection[str],
) -> None:
    """Update obscore exposure records that miss region information from
    matching visit records.

    Parameters
    ----------
    repo : `str`
        URI to the butler repository.
    check : `bool`
        If `True` then print then number of obscore records that have region
        information missing and return.
    dry_run : `bool`
        If `True` then skip final update.
    dataproduct_type : `str`
        Dataproduct type for selected exposure records, e.g. "image".
    dataproduct_subtype : `str`
        Dataproduct subtype for selected exposure records, e.g. "lsst.raw".
    instrument : `str`, optional
        Optional name of the instrument, if specified then updates are limited
        to that instrument only.
    exposure_column : `str`
        Name of the obscore column containing exposure IDs. This column is not
        a standard obscore column, but it must be configured.
    detector_column : `str`
        Name of the obscore column containing detector IDs. This column is not
        a standard obscore column, but it must be configured.
    region_columns : `Collection` [ `str` ]
        A non-empty list of column names that contain spatial information.
        These are the columns that will be updated with new data when _all_ of
        them are ``NULL``.

    Notes
    -----
    The set of updated records is limited to those that match
    ``dataproduct_type``, ``dataproduct_subtype``, ``instrument``, and
    have all values of ``region_columns`` columns set to NULL.
    """

    # Just to make sure that we have reasonable input
    assert len(region_columns) > 0, "Need at least one region-related column"

    butler = Butler(repo, writeable=True)

    # There are no client API for updating obscore table, so we are going to
    # mess with the table contents directly. For that we need to have an access
    # to the internals of the Registry and obscore manager.

    registry = butler.registry
    assert isinstance(registry, SqlRegistry), "Registry must be SqlRegistry"

    db = registry._db
    obscore_mgr = registry._managers.obscore
    if obscore_mgr is None:
        raise ValueError(f"Repository {repo} does not have obscore table.")
    assert isinstance(obscore_mgr, ObsCoreLiveTableManager), "Expect ObsCoreLiveTableManager type"

    obscore_table = obscore_mgr.table
    assert obscore_mgr.schema.dataset_fk is not None, "Live obscore table must have foreign key"
    dataset_id_column: str = obscore_mgr.schema.dataset_fk.name

    # WHERE expressions for finding all records with missing region data.
    missing_where = sqlalchemy.and_(
        obscore_table.columns["dataproduct_type"] == dataproduct_type,
        obscore_table.columns["dataproduct_subtype"] == dataproduct_subtype,
        *(obscore_table.columns[column].is_(None) for column in region_columns),
    )
    if instrument:
        missing_where = sqlalchemy.and_(obscore_table.columns["instrument_name"] == instrument, missing_where)

    def _count_missing() -> int:
        """Return count of records with missing region data."""
        query = sqlalchemy.select(sqlalchemy.func.count()).select_from(obscore_table).where(missing_where)
        with db.query(query) as result:
            return result.scalar_one()

    if check:
        # Just print count and stop here.
        count = _count_missing()
        print(
            f"Found {count} records with missing region info for" f" {dataproduct_type}/{dataproduct_subtype}"
        )
        return

    # Select all exposures with missing regions with all detectors for that
    # exposure.
    columns = (
        obscore_table.columns[dataset_id_column],
        obscore_table.columns["instrument_name"],
        obscore_table.columns[exposure_column],
        obscore_table.columns[detector_column],
    )
    # Make a list of (DatasetId, instrument, exposure, detector) for matching
    # records.
    missing_records: list[tuple[DatasetId, str, int, int]]
    query = sqlalchemy.select(*columns).select_from(obscore_table).where(missing_where)
    with db.query(query) as result:
        missing_records = [tuple(row) for row in result]  # type: ignore

    instrument_exposures: MutableMapping[str, set[int]] = defaultdict(set)
    for _, instrument, exposure, _ in missing_records:
        instrument_exposures[instrument].add(exposure)
    n_unique_exposures = sum(len(exp_set) for exp_set in instrument_exposures.values())
    _LOG.info(
        "Found %d records (%d unique exposures) with missing regions",
        len(missing_records),
        n_unique_exposures,
    )
    if not missing_records:
        return

    # Build a mapping from instrument+exposure to visit.
    exposure_to_visit: dict[tuple[str, int], int] = {}
    instrument_visits: MutableMapping[str, set[int]] = defaultdict(set)
    for instrument, exposures in instrument_exposures.items():
        # Find visit record, there may be many visits per exposure, but
        # they all should define the same region, we take any one of them.
        for exp_chunk in chunk_iterable(sorted(exposures)):
            records = registry.queryDimensionRecords(
                "visit_definition",
                where="exposure IN (exposures_param)",
                instrument=instrument,
                bind={"exposures_param": tuple(exp_chunk)},
            )
            for rec in records:
                _LOG.debug("visit_definition record=%s", rec)
                exposure_to_visit[(instrument, rec.exposure)] = rec.visit
                instrument_visits[instrument].add(rec.visit)

    # Build a mapping from a instrument+visit+detector to region
    visit_detector_region: dict[tuple[str, int, int], Region] = {}
    for instrument, visits in instrument_visits.items():
        for visit_chunk in chunk_iterable(sorted(visits)):
            records = registry.queryDimensionRecords(
                "visit_detector_region",
                where="visit IN (visits_param)",
                instrument=instrument,
                bind={"visits_param": tuple(visit_chunk)},
            )
            for rec in records:
                _LOG.debug("visit_detector_region record=%s", rec)
                visit_detector_region[(instrument, rec.visit, rec.detector)] = rec.region

    # Build dictionaries for update() method. DatasetID is the primary key in
    # obscore table, so it is sufficient to identify a record.
    update_rows: list[dict[str, Any]] = []
    for dataset_id, instrument, exposure, detector in missing_records:
        # Find a region for every visit+detector
        visit = exposure_to_visit.get((instrument, exposure))
        if visit is not None:
            region = visit_detector_region.get((instrument, visit, detector))
            if region is not None:
                region_data: Record = {}

                # ask each plugin for its values to add to a record.
                for plugin in obscore_mgr.spatial_plugins:
                    plugin_record = plugin.make_records(dataset_id, region)
                    if plugin_record is not None:
                        region_data.update(plugin_record)

                region_data["dataset_id_column"] = dataset_id
                update_rows.append(region_data)
                _LOG.debug("exposure=%s detector=%s region=%s", exposure, detector, region)

    _LOG.info("%d exposure records will be updated", len(update_rows))
    if not dry_run and update_rows:
        # Run update query
        where_dict: dict[str, Any] = {
            dataset_id_column: "dataset_id_column",
        }

        count = 0
        for rows_chunk in chunk_iterable(update_rows):
            count += db.update(obscore_table, where_dict, *rows_chunk)
        _LOG.info("Updated %s obscore records.", count)

        # Re-check number of records remaining.
        count = _count_missing()
        _LOG.info(f"{count} records with missing region info remain after update")
