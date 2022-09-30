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
from collections.abc import Collection
from typing import Any, Dict, List, Optional, Tuple

import sqlalchemy
from lsst.daf.butler import Butler, DatasetId
from lsst.daf.butler.registries.sql import SqlRegistry
from lsst.daf.butler.registry.obscore import ObsCoreLiveTableManager, RecordFactory

_LOG = logging.getLogger(__name__)


def obscore_set_exposure_regions(
    repo: str,
    check: bool,
    dry_run: bool,
    dataproduct_type: str,
    dataproduct_subtype: str,
    instrument: Optional[str],
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
        result = db.query(query)
        return result.scalar()

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
    query = sqlalchemy.select(*columns).select_from(obscore_table).where(missing_where)
    result = db.query(query)

    # Make a list of (DatasetId, instrument, exposure, detector) for matching
    # records
    missing_records: List[Tuple[DatasetId, str, int, int]] = [tuple(row) for row in result]  # type: ignore
    exposures = set((row[1], row[2]) for row in missing_records)
    _LOG.info(
        "Found %d records (%d unique exposures) with missing regions", len(missing_records), len(exposures)
    )
    if not missing_records:
        return

    # Build a mapping from exposure to visit
    exposure_to_visit: Dict[Tuple[str, int], int] = {}
    for instrument, exposure in exposures:

        # Find visit record, there may be many visits per exposure, but
        # they all should define the same region, we take any one of them.
        records = registry.queryDimensionRecords("visit_definition", instrument=instrument, exposure=exposure)
        record = next(iter(records), None)
        if record is not None:
            exposure_to_visit[(instrument, exposure)] = record.visit

    update_rows: List[Dict[str, Any]] = []
    for dataset_id, instrument, exposure, detector in missing_records:
        # Find a region for every visit+detector
        visit = exposure_to_visit.get((instrument, exposure))
        if visit is not None:
            records = registry.queryDimensionRecords(
                "visit_detector_region", instrument=instrument, visit=visit, detector=detector
            )
            record = next(iter(records), None)
            if record is not None:
                region_data: Dict[str, Any] = {}
                RecordFactory.region_to_columns(record.region, region_data)

                # Only use those columns that were given explicitly (and
                # crash if those columns are not valid).
                row: Dict[str, Any] = {}
                for column in region_columns:
                    row[column] = region_data[column]
                row["dataset_id_column"] = dataset_id
                update_rows.append(row)
                _LOG.debug("exposure=%s detector=%s region=%s", exposure, detector, record.region)

    # Build dictionaries for update() method. DatasetID is the primary key
    # in obscore table, so it is sufficient to identify a record.
    where_dict: Dict[str, Any] = {
        dataset_id_column: "dataset_id_column",
    }

    # Run update query
    if not dry_run:
        count = db.update(obscore_table, where_dict, *update_rows)
        _LOG.info("Updated %s obscore records for instrument %s", count, instrument)

        # Re-check number of records remaining.
        count = _count_missing()
        _LOG.info(f"{count} records with missing region info remain after update")
