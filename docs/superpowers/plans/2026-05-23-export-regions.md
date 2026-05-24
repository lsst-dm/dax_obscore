# `butler obscore export-regions` Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a new `butler obscore export-regions` subcommand that writes a FITS binary table of distinct STC-S region strings for the datasets matched by a collection / where / dataset-type query, suitable as input to TOPCAT STILTS `mocshape`.

**Architecture:** A thin click subcommand inside the existing `obscore` group delegates to a script function `obscore_export_regions(...)`. The script opens the butler, resolves the region dimension via `DimensionGroup.region_dimension`, runs a `butler.query().join_dataset_search().where().general()` query that yields the dataId and region per matching dataset, deduplicates by the region dimension's required keys (sorted for stability), converts each region to STC-S via `lsst.sphgeom.Region.to_ivoa_stcs()`, and writes a single-extension FITS BinTable with a `s_region` column tagged with `TXTYP1='stc-s'` so STILTS auto-detects it.

**Tech Stack:** Python 3, click, `lsst.daf.butler`, `lsst.sphgeom` (`tickets/DM-53569` branch), `astropy.io.fits`, unittest.

**Spec:** `docs/superpowers/specs/2026-05-23-export-regions-design.md`

---

## File Structure

**Created:**

- `python/lsst/dax/obscore/script/obscore_export_regions.py` — `obscore_export_regions(repo, destination, dataset_type, collections, where)` script function.
- `tests/test_export_regions.py` — unit tests for the script function plus a CLI smoke test.
- `doc/changes/DM-53569.feature` — release note fragment.

**Modified:**

- `python/lsst/dax/obscore/script/__init__.py` — re-export `obscore_export_regions`.
- `python/lsst/dax/obscore/cli/cmd/commands.py` — register the new `export-regions` click subcommand inside the `obscore` group.

---

## Task 1: Script module skeleton

**Files:**
- Create: `python/lsst/dax/obscore/script/obscore_export_regions.py`
- Modify: `python/lsst/dax/obscore/script/__init__.py`

Goal: a minimal, importable function so subsequent test tasks can import it. The body raises `NotImplementedError` so that the first failing test fails for the *right* reason.

- [ ] **Step 1: Create the script module with a stub function**

Create `python/lsst/dax/obscore/script/obscore_export_regions.py`:

```python
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

from collections.abc import Iterable


def obscore_export_regions(
    repo: str,
    destination: str,
    dataset_type: str,
    collections: Iterable[str],
    where: str,
) -> None:
    """Export distinct spatial regions for the datasets matching a query
    as a FITS binary table of IVOA STC-S strings.

    Parameters
    ----------
    repo : `str`
        URI to the butler repository.
    destination : `str`
        Path to the output FITS file.
    dataset_type : `str`
        Name of the butler dataset type whose regions should be exported.
    collections : `~collections.abc.Iterable` [ `str` ]
        Collections to search.
    where : `str`
        Butler query ``where`` expression. Empty string means no constraint.
    """
    raise NotImplementedError
```

- [ ] **Step 2: Re-export the function from `script/__init__.py`**

Add an import line to `python/lsst/dax/obscore/script/__init__.py` next to the existing `obscore_*` imports:

```python
from .obscore_export_regions import *
```

(The existing file uses `from .obscore_export import *` style — match that.)

- [ ] **Step 3: Verify the import works**

Run: `python -c "from lsst.dax.obscore.script import obscore_export_regions; print(obscore_export_regions)"`

Expected: prints the function reference, no exception.

- [ ] **Step 4: Commit**

```bash
git add python/lsst/dax/obscore/script/obscore_export_regions.py \
        python/lsst/dax/obscore/script/__init__.py
git commit -m "Add obscore_export_regions stub for DM-53569"
```

---

## Task 2: Test scaffolding + first failing test (visit/detector happy path)

**Files:**
- Create: `tests/test_export_regions.py`

Reuses the existing test fixture (`DaxObsCoreTestMixin.make_butler` + `tests/data/hsc_gen3.yaml`) which loads both `_mock_calexp` (visit + detector dataset type) and `_mock_deepCoadd` (skymap + tract + patch dataset type).

- [ ] **Step 1: Write the failing test**

Create `tests/test_export_regions.py`:

```python
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

import os
import unittest

import astropy.io.fits as fits

from lsst.daf.butler.tests.utils import makeTestTempDir, removeTestTempDir
from lsst.dax.obscore.script import obscore_export_regions
from lsst.dax.obscore.tests import DaxObsCoreTestMixin

TESTDIR = os.path.abspath(os.path.dirname(__file__))


class ExportRegionsTestCase(unittest.TestCase, DaxObsCoreTestMixin):
    """Tests for obscore_export_regions."""

    def setUp(self):
        self.root = makeTestTempDir(TESTDIR)

    def tearDown(self):
        removeTestTempDir(self.root)

    def _populate(self):
        butler = self.make_butler()
        self.enterContext(butler)
        butler.import_(
            filename=os.path.join(TESTDIR, "data", "hsc_gen3.yaml"),
            without_datastore=True,
        )
        return butler

    def test_visit_detector_happy_path(self):
        """Distinct visit_detector_region regions are written as STC-S."""
        butler = self._populate()
        out = os.path.join(self.root, "regions.fits")

        obscore_export_regions(
            repo=self.root,
            destination=out,
            dataset_type="_mock_calexp",
            collections=["HSC/runs/ci_hsc"],
            where="",
        )

        self.assertTrue(os.path.exists(out))
        with fits.open(out) as hdul:
            self.assertEqual(len(hdul), 2)  # primary + bintable
            tbl = hdul[1]
            self.assertEqual(tbl.header["EXTNAME"], "REGIONS")
            self.assertEqual(tbl.columns.names, ["s_region"])
            self.assertEqual(tbl.header["TXTYP1"], "stc-s")
            self.assertEqual(
                tbl.header["TUTYP1"],
                "obscore:Char.SpatialAxis.Coverage.Support.Area",
            )
            self.assertEqual(tbl.header["TUCD1"], "pos.outline;obs.field")
            data = tbl.data["s_region"]
            self.assertGreater(len(data), 0)
            for s in data:
                # Every value parses as STC-S (one of the allowed leading
                # tokens). Frame keyword is always second.
                head = s.split()[0]
                self.assertIn(
                    head,
                    {"Circle", "Polygon", "Ellipse", "Union", "Intersection"},
                )
                self.assertEqual(s.split()[1], "ICRS")


if __name__ == "__main__":
    unittest.main()
```

Notes on choices:

- `make_butler()` creates the repo at `self.root`, so `repo=self.root` is the simplest URI to pass.
- The fixture's collection name is `HSC/runs/ci_hsc` (matches `make_export_config` and `hsc_gen3.yaml`).

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_export_regions.py::ExportRegionsTestCase::test_visit_detector_happy_path -v`

Expected: FAIL with `NotImplementedError` raised from the stub.

- [ ] **Step 3: Commit the failing test**

```bash
git add tests/test_export_regions.py
git commit -m "Add failing test for obscore_export_regions"
```

---

## Task 3: Implement the script function

Now implement the real `obscore_export_regions` so the test from Task 2 passes. This is the largest implementation task; subsequent tests mostly verify edge cases the implementation already handles.

**Files:**
- Modify: `python/lsst/dax/obscore/script/obscore_export_regions.py`

- [ ] **Step 1: Replace the stub with the full implementation**

Replace the file contents with:

```python
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
from collections.abc import Iterable

import astropy.io.fits as fits
import numpy as np

from lsst.daf.butler import Butler
from lsst.sphgeom import Region

_LOG = logging.getLogger(__name__)


def obscore_export_regions(
    repo: str,
    destination: str,
    dataset_type: str,
    collections: Iterable[str],
    where: str,
) -> None:
    """Export distinct spatial regions for the datasets matching a query
    as a FITS binary table of IVOA STC-S strings.

    Parameters
    ----------
    repo : `str`
        URI to the butler repository.
    destination : `str`
        Path to the output FITS file.
    dataset_type : `str`
        Name of the butler dataset type whose regions should be exported.
    collections : `~collections.abc.Iterable` [ `str` ]
        Collections to search.
    where : `str`
        Butler query ``where`` expression. Empty string means no constraint.
    """
    collections = list(collections)
    with Butler.from_config(repo, writeable=False) as butler:
        dt = butler.get_dataset_type(dataset_type)
        region_dim = dt.dimensions.region_dimension
        if region_dim is None:
            raise ValueError(
                f"Dataset type {dataset_type!r} has no spatial region dimension."
            )
        region_key = f"{region_dim}.region"
        required_keys = sorted(
            butler.dimensions[region_dim].minimal_group.required
        )

        regions: list[Region] = []
        seen: set[tuple] = set()
        with butler.query() as query:
            q = query.join_dataset_search(dt.name, collections=collections)
            if where:
                q = q.where(where)
            result = q.general(
                dt.dimensions,
                region_key,
                dataset_fields={dt.name: ...},
                find_first=True,
            )
            for dataId, _refs, raw_row in result.iter_tuples(dt):
                key = tuple(dataId[k] for k in required_keys)
                if key in seen:
                    continue
                seen.add(key)
                regions.append(raw_row[region_key])

    stcs_strings = [region.to_ivoa_stcs() for region in regions]
    _write_fits(stcs_strings, destination, dataset_type, collections, where)


def _write_fits(
    stcs_strings: list[str],
    destination: str,
    dataset_type: str,
    collections: list[str],
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
        f"export-regions dataset_type={dataset_type} "
        f"collections={','.join(collections)} where={where!r}"
    )
    hdul = fits.HDUList([fits.PrimaryHDU(), hdu])
    hdul.writeto(destination, overwrite=True)
```

Notes:

- Fixed-width `An` column (rather than variable-length `PA`) is used because astropy's variable-length string handling is awkward and STILTS reads `An` cleanly. Per the spec, either is acceptable.
- `HISTORY` line records the query for traceability — same idea as the SIAv2 `query_url` field.

- [ ] **Step 2: Run the Task 2 test to verify it passes**

Run: `pytest tests/test_export_regions.py::ExportRegionsTestCase::test_visit_detector_happy_path -v`

Expected: PASS.

- [ ] **Step 3: Commit**

```bash
git add python/lsst/dax/obscore/script/obscore_export_regions.py
git commit -m "Implement obscore_export_regions"
```

---

## Task 4: Test patch (skymap+tract+patch) happy path

**Files:**
- Modify: `tests/test_export_regions.py`

- [ ] **Step 1: Add the test method**

Append the following method to `ExportRegionsTestCase`:

```python
    def test_patch_happy_path(self):
        """Distinct patch regions are written as STC-S."""
        butler = self._populate()
        out = os.path.join(self.root, "regions.fits")

        obscore_export_regions(
            repo=self.root,
            destination=out,
            dataset_type="_mock_deepCoadd",
            collections=["HSC/runs/ci_hsc"],
            where="",
        )

        with fits.open(out) as hdul:
            data = hdul[1].data["s_region"]
            self.assertGreater(len(data), 0)
            for s in data:
                head = s.split()[0]
                self.assertIn(
                    head,
                    {"Circle", "Polygon", "Ellipse", "Union", "Intersection"},
                )
```

- [ ] **Step 2: Run the test**

Run: `pytest tests/test_export_regions.py::ExportRegionsTestCase::test_patch_happy_path -v`

Expected: PASS (the implementation already supports this case).

- [ ] **Step 3: Commit**

```bash
git add tests/test_export_regions.py
git commit -m "Test patch dataset type for export-regions"
```

---

## Task 5: Test deduplication

Confirms the dedup logic actually deduplicates — i.e., row count matches the number of distinct `(instrument, visit, detector)` tuples in the matched datasets, not the number of datasets.

**Files:**
- Modify: `tests/test_export_regions.py`

- [ ] **Step 1: Add the test method**

```python
    def test_dedup_visit_detector(self):
        """Row count equals distinct (instrument, visit, detector) keys."""
        butler = self._populate()
        out = os.path.join(self.root, "regions.fits")

        # Build the expected distinct-key count from the dataset list
        # itself, using the dataId of every matching ref. This avoids
        # depending on butler query APIs that may not be relevant here.
        refs = list(
            butler.query_datasets(
                "_mock_calexp", collections=["HSC/runs/ci_hsc"]
            )
        )
        expected_keys = {
            (ref.dataId["instrument"], ref.dataId["visit"], ref.dataId["detector"])
            for ref in refs
        }
        # Sanity check: the fixture must have duplicates for this test
        # to actually exercise dedup.
        self.assertGreater(
            len(refs),
            len(expected_keys),
            "Fixture must have multiple datasets per (visit, detector) "
            "for the dedup test to be meaningful",
        )

        obscore_export_regions(
            repo=self.root,
            destination=out,
            dataset_type="_mock_calexp",
            collections=["HSC/runs/ci_hsc"],
            where="",
        )

        with fits.open(out) as hdul:
            self.assertEqual(len(hdul[1].data), len(expected_keys))
```

Notes:

- `butler.query_datasets` is the modern API used elsewhere in this repo (commit `59ad89d`).
- If the fixture happens not to have duplicates, the `assertGreater` will fail and tell you to pick a different dataset type or amend the fixture. We do not silently skip in that case.

- [ ] **Step 2: Run the test**

Run: `pytest tests/test_export_regions.py::ExportRegionsTestCase::test_dedup_visit_detector -v`

Expected: PASS.

- [ ] **Step 3: Commit**

```bash
git add tests/test_export_regions.py
git commit -m "Test region deduplication for export-regions"
```

---

## Task 6: Test `--where` filtering

**Files:**
- Modify: `tests/test_export_regions.py`

- [ ] **Step 1: Add the test method**

```python
    def test_where_filtering(self):
        """A where clause restricts the produced regions."""
        butler = self._populate()
        out_all = os.path.join(self.root, "regions_all.fits")
        out_one = os.path.join(self.root, "regions_one.fits")

        # Pick one (instrument, visit, detector) triple from the fixture
        # by enumerating the matching dataset refs.
        refs = list(
            butler.query_datasets(
                "_mock_calexp", collections=["HSC/runs/ci_hsc"]
            )
        )
        keys = sorted(
            {
                (ref.dataId["instrument"], ref.dataId["visit"], ref.dataId["detector"])
                for ref in refs
            }
        )
        self.assertGreater(
            len(keys), 1,
            "Need >1 distinct (visit, detector) pair in the fixture",
        )
        instrument, visit, detector = keys[0]
        where = (
            f"instrument='{instrument}' "
            f"AND visit={visit} "
            f"AND detector={detector}"
        )

        obscore_export_regions(
            repo=self.root,
            destination=out_all,
            dataset_type="_mock_calexp",
            collections=["HSC/runs/ci_hsc"],
            where="",
        )
        obscore_export_regions(
            repo=self.root,
            destination=out_one,
            dataset_type="_mock_calexp",
            collections=["HSC/runs/ci_hsc"],
            where=where,
        )

        with fits.open(out_all) as hdul_all, fits.open(out_one) as hdul_one:
            self.assertEqual(len(hdul_one[1].data), 1)
            self.assertLess(len(hdul_one[1].data), len(hdul_all[1].data))
```

- [ ] **Step 2: Run the test**

Run: `pytest tests/test_export_regions.py::ExportRegionsTestCase::test_where_filtering -v`

Expected: PASS.

- [ ] **Step 3: Commit**

```bash
git add tests/test_export_regions.py
git commit -m "Test where filtering for export-regions"
```

---

## Task 7: Test error when dataset type has no spatial dimension

**Files:**
- Modify: `tests/test_export_regions.py`

The test fixture should expose at least one dataset type with no spatial dimension (e.g., a calibration-style type). If none is convenient, we register one inline using the butler's API.

- [ ] **Step 1: Add the test method**

```python
    def test_error_no_region_dimension(self):
        """A dataset type with no region dimension is rejected."""
        from lsst.daf.butler import DatasetType

        butler = self._populate()
        # Register a dataset type whose dimensions have no region.
        # `instrument` alone has no region dimension.
        dt = DatasetType(
            "_mock_no_region",
            dimensions=("instrument",),
            storageClass="StructuredDataDict",
            universe=butler.dimensions,
        )
        butler.registry.registerDatasetType(dt)

        out = os.path.join(self.root, "regions.fits")
        with self.assertRaisesRegex(ValueError, "no spatial region dimension"):
            obscore_export_regions(
                repo=self.root,
                destination=out,
                dataset_type="_mock_no_region",
                collections=["HSC/runs/ci_hsc"],
                where="",
            )
```

Notes:

- `StructuredDataDict` is a stable storage class name in daf_butler that requires no extra setup. If it is not available in the test universe, substitute any storage class that exists in the test repo's storage class registry — e.g., `"DataFrame"` or look at `butler.storageClasses.keys()` once during local development.

- [ ] **Step 2: Run the test**

Run: `pytest tests/test_export_regions.py::ExportRegionsTestCase::test_error_no_region_dimension -v`

Expected: PASS.

- [ ] **Step 3: Commit**

```bash
git add tests/test_export_regions.py
git commit -m "Test error for non-spatial dataset type in export-regions"
```

---

## Task 8: Test empty result

**Files:**
- Modify: `tests/test_export_regions.py`

Confirms an empty result still produces a valid FITS file (zero rows), not an exception or a malformed file.

- [ ] **Step 1: Add the test method**

```python
    def test_empty_result(self):
        """A where clause that matches nothing still produces a valid FITS file."""
        butler = self._populate()  # noqa: F841 -- ensures repo populated
        out = os.path.join(self.root, "regions_empty.fits")

        obscore_export_regions(
            repo=self.root,
            destination=out,
            dataset_type="_mock_calexp",
            collections=["HSC/runs/ci_hsc"],
            where="visit=-1",
        )

        with fits.open(out) as hdul:
            self.assertEqual(len(hdul[1].data), 0)
            self.assertEqual(hdul[1].columns.names, ["s_region"])
            self.assertEqual(hdul[1].header["TXTYP1"], "stc-s")
```

- [ ] **Step 2: Run the test**

Run: `pytest tests/test_export_regions.py::ExportRegionsTestCase::test_empty_result -v`

Expected: PASS.

- [ ] **Step 3: Commit**

```bash
git add tests/test_export_regions.py
git commit -m "Test empty result handling for export-regions"
```

---

## Task 9: Wire up the click subcommand

**Files:**
- Modify: `python/lsst/dax/obscore/cli/cmd/commands.py`

- [ ] **Step 1: Add the new subcommand to `commands.py`**

Append a new `@obscore.command` block at the end of `python/lsst/dax/obscore/cli/cmd/commands.py`. Place it after `def siav2(...)`.

```python
@obscore.command(
    short_help="Export distinct spatial regions for matching datasets as STC-S.",
    cls=ButlerCommand,
)
@repo_argument(required=True)
@destination_argument(
    required=True,
    help="DESTINATION is the location of the output FITS file.",
    type=MWPath(file_okay=True, dir_okay=False, writable=True),
)
@click.option(
    "--dataset-type",
    "dataset_type",
    help="Name of the butler dataset type whose regions should be exported.",
    required=True,
    type=str,
)
@collections_option()
@where_option()
@options_file_option()
def export_regions(*args: Any, **kwargs: Any) -> None:
    """Export distinct spatial regions for the datasets matching a query
    as a FITS binary table of IVOA STC-S strings.

    The output is suitable as input to the TOPCAT STILTS ``mocshape``
    command. The dataset type's ``DimensionGroup.region_dimension`` is used
    to decide which spatial dimension supplies the region, naturally
    preferring ``visit_detector_region`` over ``visit`` and ``patch`` over
    ``tract``.
    """
    script.obscore_export_regions(*args, **kwargs)
```

Notes:

- `dataset_type` is declared inline (single string, required) rather than via the comma-splitting `dataset_type_option` shared by `export` and `siav2`. The plural form is wrong here — we want one dataset type only.
- The Click parameter is exposed as `dataset_type` (snake-case) so it maps directly onto the script function's parameter name.

- [ ] **Step 2: Verify the command appears in the CLI**

Run: `butler obscore --help`

Expected: the output lists `export-regions` alongside `export`, `set-exposure-regions`, `update-table`, `siav2`.

- [ ] **Step 3: Verify `--help` for the new subcommand**

Run: `butler obscore export-regions --help`

Expected: shows `--dataset-type`, `--collections`, `--where`, repo and destination arguments. Exit code 0.

- [ ] **Step 4: Commit**

```bash
git add python/lsst/dax/obscore/cli/cmd/commands.py
git commit -m "Register obscore export-regions click subcommand"
```

---

## Task 10: CLI smoke test

**Files:**
- Modify: `tests/test_export_regions.py`

End-to-end test through the Click CLI to confirm the wiring works.

- [ ] **Step 1: Add the smoke test**

Append the following to `tests/test_export_regions.py`. The import goes near the top with the other imports:

```python
from click.testing import CliRunner

from lsst.dax.obscore.cli.cmd.commands import obscore as obscore_cli
```

And the test method on `ExportRegionsTestCase`:

```python
    def test_cli_smoke(self):
        """The click subcommand wiring works end-to-end."""
        self._populate()
        out = os.path.join(self.root, "regions_cli.fits")
        runner = CliRunner()
        result = runner.invoke(
            obscore_cli,
            [
                "export-regions",
                self.root,
                out,
                "--dataset-type",
                "_mock_calexp",
                "--collections",
                "HSC/runs/ci_hsc",
            ],
        )
        if result.exit_code != 0:
            # Surface the exception for easier debugging.
            raise AssertionError(
                f"exit_code={result.exit_code}\n"
                f"output={result.output}\n"
                f"exception={result.exception!r}"
            )
        self.assertTrue(os.path.exists(out))
        with fits.open(out) as hdul:
            self.assertGreater(len(hdul[1].data), 0)
```

- [ ] **Step 2: Run the test**

Run: `pytest tests/test_export_regions.py::ExportRegionsTestCase::test_cli_smoke -v`

Expected: PASS.

- [ ] **Step 3: Run the entire test module**

Run: `pytest tests/test_export_regions.py -v`

Expected: all tests in the module pass.

- [ ] **Step 4: Commit**

```bash
git add tests/test_export_regions.py
git commit -m "Add CLI smoke test for export-regions"
```

---

## Task 11: Release note fragment

**Files:**
- Create: `doc/changes/DM-53569.feature`

- [ ] **Step 1: Write the news fragment**

Create `doc/changes/DM-53569.feature`:

```
Added ``butler obscore export-regions`` subcommand that writes a FITS binary table of distinct IVOA STC-S region strings for the datasets matching a collection / where / dataset-type query. Suitable as input to the TOPCAT STILTS ``mocshape`` command.
```

- [ ] **Step 2: Commit**

```bash
git add doc/changes/DM-53569.feature
git commit -m "Add news fragment for DM-53569"
```

---

## Task 12: Final verification

- [ ] **Step 1: Run the full test suite**

Run: `pytest tests/ -v`

Expected: all tests pass. No new failures introduced in other modules.

- [ ] **Step 2: Run mypy and pre-commit hooks**

Run: `pre-commit run --all-files`

Expected: clean. If any check fails, fix the underlying issue and amend (or follow up with a fixup commit) — do not skip hooks.

- [ ] **Step 3: Manual sanity check (optional but recommended)**

If a real butler repo is available, run:

```
butler obscore export-regions <repo> /tmp/regions.fits \
    --dataset-type calexp --collections <collection> --where "instrument='LSSTCam'"
```

Then:

```
stilts mocshape in=/tmp/regions.fits coords=s_region order=10 out=/tmp/coverage.moc.fits
```

Expected: STILTS reports it auto-detected `shape=stc-s` and writes a MOC file.
