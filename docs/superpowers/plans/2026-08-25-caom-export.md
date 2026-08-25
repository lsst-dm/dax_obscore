# CAOM Export Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a `butler obscore export-caom` command that writes CAOM-2.4 XML documents from a Butler repository, driven by a `caom` section in the existing ObsCore configuration files.

**Architecture:** `CaomExporter` consumes the same ObsCore records that `ObscoreExporter` already produces, so the CAOM output cannot drift from the ObsCore output. One contained refactor extracts the query loop from `ObscoreExporter._make_record_batches` into a shared generator. Records are grouped into `Observation` → `Plane` → `Artifact` using per-dataset-type templates, and written one XML document per Observation. Nothing reads a data file.

**Tech Stack:** Python 3.11+, pydantic v2, `lsst.daf.butler`, `lsst.sphgeom`, `astropy.units`, `caom2` (optional dependency), click, pytest.

**Spec:** `docs/superpowers/specs/2026-08-25-caom-export-design.md`

## Global Constraints

- Python 3.11+. `from __future__ import annotations` at the top of every new module.
- Every new file starts with the GPL header block copied verbatim from `python/lsst/dax/obscore/config.py` lines 1-20, with "This file is part of dax_obscore." as the first line.
- ruff and mypy must both pass. `line-length = 110`, `max-doc-length = 79`, numpydoc docstring convention, `known-first-party = ["lsst"]`.
- numpydoc validation runs in pre-commit. Every public function needs Parameters and Returns sections.
- **All unit conversion goes through `astropy.units`. Do not write `3600`, `math.pi/180`, or any other conversion constant anywhere in this work.**
- `caom2` is an optional dependency. `caom_config.py` must never import it, directly or transitively. Only `caom_exporter.py` and `_caom_shapes.py` may, behind the lazy guard from Task 2.
- Never reformat or restructure code outside the lines a task names.
- Run `pytest tests/ -x -q` before every commit. All existing tests must keep passing.

## Phases

Tasks 1-9 deliver phases 1 to 3 of the spec. Phase 4 is a review round with CADC and has no code. Phase 5 is a companion `daf_butler` ticket moving `s_pixel_scale` onto `DatasetTypeConfig`; it is described in "Follow-on work" at the end and is deliberately not planned in detail, because phase 4 may change it.

- Task 1 is spec phase 1.
- Tasks 2-8 are spec phase 2.
- Tasks 9-10 are spec phase 3.

## File Structure

| File | Responsibility |
| --- | --- |
| `python/lsst/dax/obscore/obscore_exporter.py` | Modified: query loop extracted into `_iter_record_refs`; new `_resolve_uris` helper |
| `python/lsst/dax/obscore/caom_config.py` | Pydantic models for the `caom` block. No `caom2` import, ever |
| `python/lsst/dax/obscore/_caom_shapes.py` | Pure conversion helpers: sphgeom region to caom2 shape, unit conversion |
| `python/lsst/dax/obscore/caom_exporter.py` | `CaomExporter`: grouping, CAOM object construction, XML output |
| `python/lsst/dax/obscore/script/obscore_export_caom.py` | Script layer |
| `python/lsst/dax/obscore/cli/cmd/commands.py` | Modified: new `export-caom` subcommand |
| `python/lsst/dax/obscore/config.py` | Modified: `ExporterConfig.caom` field |
| `python/lsst/dax/obscore/tests.py` | Modified: `make_caom_config` helper for tests |
| `tests/test_caom_config.py` | Config validation. Runs without `caom2` |
| `tests/test_caom_shapes.py` | Region and unit conversion. Requires `caom2` |
| `tests/test_caom_exporter.py` | End-to-end export against the mock butler. Requires `caom2` |

---

### Task 1: Extract the shared record generator

Pure refactor of `ObscoreExporter`. No output changes anywhere. This is spec phase 1 and must be committed on its own so that any regression is bisectable.

**Files:**
- Modify: `python/lsst/dax/obscore/obscore_exporter.py` (the `_make_record_batches` method, currently the last method in the file)
- Test: `tests/test_exporter.py` (existing, unmodified)

**Interfaces:**
- Consumes: nothing.
- Produces: `_QueryState` dataclass with a single `overflow: bool = False` field. `ObscoreExporter._iter_record_refs(state: _QueryState, limit: int | None = None) -> Iterator[tuple[DatasetRef, Region | None, dict[str, Any]]]`.

- [ ] **Step 1: Run the existing tests to record the baseline**

Run: `pytest tests/test_exporter.py -q`

Expected: all pass. Note the count; it must be identical after the refactor.

- [ ] **Step 2: Add the `_QueryState` dataclass**

Add `import dataclasses` to the imports at the top of `obscore_exporter.py`, keeping isort order (it goes with the other stdlib imports, after `contextlib`).

Add this class immediately before `class ObscoreExporter:`:

```python
@dataclasses.dataclass
class _QueryState:
    """Mutable state shared between the record generator and its callers.

    Notes
    -----
    The generator cannot return a value to a caller that is iterating it,
    so the overflow flag is carried here instead.
    """

    overflow: bool = False
    """`True` if more records matched the query than the requested limit."""
```

- [ ] **Step 3: Add `_iter_record_refs`**

Add this method to `ObscoreExporter`, immediately before the existing `_make_record_batches`. The body is the existing `_make_record_batches` body with the batching removed, `overflow` replaced by `state.overflow`, and `yield ref, region, record` in place of `batch.add_to_batch(record)`.

```python
    def _iter_record_refs(
        self, state: _QueryState, limit: int | None = None
    ) -> Iterator[tuple[DatasetRef, Region | None, dict[str, Any]]]:
        """Query the registry and generate ObsCore records with their refs.

        Parameters
        ----------
        state : `_QueryState`
            Mutable state updated in place. The ``overflow`` attribute is
            set to `True` if more records matched than ``limit`` allowed.
        limit : `int` or `None`, optional
            Maximum number of records to generate. If `None` there is no
            limit.

        Yields
        ------
        ref : `~lsst.daf.butler.DatasetRef`
            Reference to the dataset the record describes.
        region : `~lsst.sphgeom.Region` or `None`
            Spatial region associated with the dataset, if the dataset type
            has a relevant spatial dimension.
        record : `dict` [ `str`, `~typing.Any` ]
            The ObsCore record, keyed by column name.
        """
        collections: Any = self.config.collections
        if not collections:
            raise ValueError("No collections specified. Querying all collections is not allowed.")

        if limit is not None:
            if limit == 0:
                # Return immediately since no records requested.
                return
            # Always ask for one extra to allow overflow detection.
            limit = abs(limit) + 1

        for dataset_type_name in self.config.dataset_types:
            _LOG.verbose("Querying datasets for dataset type %s [limit=%s]", dataset_type_name, limit)
            where_clauses = self.config.dataset_type_constraints.get(dataset_type_name, [self.config.where])
            if not where_clauses:
                # Want an empty default to match everything.
                where_clauses = [WhereBind(where="")]

            # Determine the relevant dimension for the region that can be
            # joined by the query system.
            dataset_type = self.butler.get_dataset_type(dataset_type_name)
            region_dim, region_metadata_name = self.record_factory.region_dimension(dataset_type.dimensions)
            region_key: str | None = None
            if region_dim is not None:
                region_key = f"{region_dim}.{region_metadata_name}"

            with self.butler.query() as query:
                for where_clause in where_clauses:
                    where_query = query

                    if where_clause.extra_dims:
                        where_query = where_query.join_dimensions(where_clause.extra_dims)

                    where_query = where_query.join_dataset_search(dataset_type_name, collections=collections)

                    if where_clause.where:
                        _LOG.verbose("Processing query with constraint %s", where_clause)
                        where_query = where_query.where(where_clause.where, bind=where_clause.bind)

                    region_args = [region_key] if region_key else []
                    result = where_query.general(
                        dataset_type.dimensions,
                        *region_args,
                        dataset_fields={dataset_type_name: ...},
                        find_first=True,
                    )

                    # We need dimension records.
                    result = result.with_dimension_records()

                    if limit is not None:
                        result = result.limit(limit)

                    count = 0
                    for dataId, (ref,), raw_row in result.iter_tuples(dataset_type):
                        dataId = ref.dataId
                        region = raw_row[region_key] if region_key else None
                        _LOG.debug("New record, dataId=%s region=%s", dataId.mapping, region)

                        self._derived_region_factory.set(dataId, region)
                        record = self.record_factory(ref)
                        if record is None:
                            continue

                        count += 1
                        if limit is not None and count == limit:
                            # Hit the +1 so should not add this to the batch.
                            _LOG.debug("Got one more than requested limit so dropping final record.")
                            state.overflow = True
                            break

                        yield ref, region, record

                    if limit is not None:
                        limit -= count
                    if state.overflow:
                        # We counted one too many so adjust for the log
                        # message.
                        count -= 1

                    _LOG.info("Copied %d records from dataset type %s", count, dataset_type_name)

                    if state.overflow:
                        # No more queries need to run.
                        # This breaks out one level of nesting.
                        break

                if state.overflow:
                    # Stop further dataset type queries.
                    break
```

- [ ] **Step 4: Replace the body of `_make_record_batches`**

Replace the entire existing `_make_record_batches` body, keeping its signature and docstring, with:

```python
    def _make_record_batches(
        self, batch_size: int = 10_000, limit: int | None = None
    ) -> Iterator[tuple[RecordBatch, bool]]:
        """Generate batches of records to save to a file.

        Yields the batches and a flag indicating whether an overflow condition
        was hit.
        """
        batch = _BatchCollector(self.schema)
        state = _QueryState()

        for _ref, _region, record in self._iter_record_refs(state, limit=limit):
            batch.add_to_batch(record)
            if batch.size >= batch_size:
                _LOG.debug("Saving next record batch, size=%s", batch.size)
                yield (batch.make_record_batch(), state.overflow)

        # Final batch if anything is there
        if batch.size > 0:
            _LOG.debug("Saving final record batch, size=%s", batch.size)
            yield (batch.make_record_batch(), state.overflow)
```

Add `DatasetRef` to the existing `from lsst.daf.butler import Butler, DataCoordinate, ddl` import line, keeping alphabetical order: `from lsst.daf.butler import Butler, DataCoordinate, DatasetRef, ddl`.

- [ ] **Step 5: Verify no behavior changed**

Run: `pytest tests/test_exporter.py tests/test_siav2.py -q`

Expected: identical pass count to Step 1. `test_export_parquet` must still report 31 columns and 35 rows, and the SIAv2 overflow tests must still pass — those exercise the `state.overflow` path.

- [ ] **Step 6: Lint**

Run: `ruff check python/lsst/dax/obscore/obscore_exporter.py && ruff format --check python/lsst/dax/obscore/obscore_exporter.py && mypy python/lsst/dax/obscore/obscore_exporter.py`

Expected: clean.

- [ ] **Step 7: Commit**

```bash
git add python/lsst/dax/obscore/obscore_exporter.py
git commit -m "Extract shared record generator from ObscoreExporter

_make_record_batches becomes a batching wrapper over _iter_record_refs,
which yields (ref, region, record) so other exporters can consume the
same query. No change to any existing output."
```

---

### Task 2: Optional `caom2` dependency and lazy import guard

Establishes the packaging and the import guard before any CAOM code exists, so every later task has somewhere safe to import from.

**Files:**
- Create: `python/lsst/dax/obscore/_caom_shapes.py`
- Create: `requirements-caom.txt`
- Modify: `pyproject.toml` (the `[project.optional-dependencies]` block, currently containing `postgres` and `test`)
- Modify: `.github/workflows/build.yaml`
- Test: `tests/test_caom_shapes.py`

**Interfaces:**
- Consumes: nothing.
- Produces: `_caom_shapes.caom2` (the module or `None`), `_caom_shapes.require_caom2() -> None` which raises `ImportError` with an actionable message when `caom2` is absent.

- [ ] **Step 1: Add the optional dependency**

In `pyproject.toml`, inside `[project.optional-dependencies]`, add a `caom` entry alongside the existing `postgres` entry:

```toml
[project.optional-dependencies]
postgres = ["psycopg2"]
caom = ["caom2 >= 2.6"]

test = [
    "pytest >= 3.2",
]
```

Create `requirements-caom.txt`:

```
caom2 >= 2.6
```

- [ ] **Step 2: Install it locally**

Run: `pip install -r requirements-caom.txt`

Expected: `caom2` installs. Confirm with `python -c "import caom2; print(caom2.__file__)"`.

- [ ] **Step 3: Write the failing test**

Create `tests/test_caom_shapes.py` with the GPL header, then:

```python
import unittest

from lsst.dax.obscore import _caom_shapes


class ImportGuardTestCase(unittest.TestCase):
    """Tests of the optional caom2 import guard."""

    def test_require_caom2_succeeds_when_available(self):
        """require_caom2 is a no-op when the library is installed."""
        self.assertIsNotNone(_caom_shapes.caom2)
        _caom_shapes.require_caom2()

    def test_require_caom2_message_names_the_extra(self):
        """The failure message tells the user how to fix the problem."""
        message = _caom_shapes._CAOM2_IMPORT_MESSAGE
        self.assertIn("lsst-dax-obscore[caom]", message)


if __name__ == "__main__":
    unittest.main()
```

- [ ] **Step 4: Run it to verify it fails**

Run: `pytest tests/test_caom_shapes.py -q`

Expected: FAIL with `ModuleNotFoundError: No module named 'lsst.dax.obscore._caom_shapes'`.

- [ ] **Step 5: Create the module with the guard**

Create `python/lsst/dax/obscore/_caom_shapes.py` with the GPL header, then:

```python
from __future__ import annotations

__all__ = ["require_caom2"]

try:
    import caom2
except ImportError:
    caom2 = None  # type: ignore[assignment]

_CAOM2_IMPORT_MESSAGE = (
    "The caom2 package is required for CAOM export but is not installed. "
    "Install it with 'pip install lsst-dax-obscore[caom]'."
)


def require_caom2() -> None:
    """Check that the optional ``caom2`` dependency is importable.

    Raises
    ------
    ImportError
        Raised if ``caom2`` is not installed.
    """
    if caom2 is None:
        raise ImportError(_CAOM2_IMPORT_MESSAGE)
```

- [ ] **Step 6: Run the test to verify it passes**

Run: `pytest tests/test_caom_shapes.py -q`

Expected: 2 passed.

- [ ] **Step 7: Install the extra in CI**

In `.github/workflows/build.yaml`, find the step that installs `requirements.txt` and add `-r requirements-caom.txt` to the same `pip install` invocation, so the CAOM tests run rather than skip.

- [ ] **Step 8: Lint and commit**

Run: `ruff check python/lsst/dax/obscore/_caom_shapes.py tests/test_caom_shapes.py && mypy python/lsst/dax/obscore/_caom_shapes.py`

```bash
git add pyproject.toml requirements-caom.txt .github/workflows/build.yaml \
    python/lsst/dax/obscore/_caom_shapes.py tests/test_caom_shapes.py
git commit -m "Add optional caom2 dependency with a lazy import guard"
```

---

### Task 3: Region and unit conversion helpers

Pure functions with no Butler and no configuration, so they are cheap to test exhaustively. This is where the spec's "no hand-written conversion constants" rule is enforced.

**Files:**
- Modify: `python/lsst/dax/obscore/_caom_shapes.py`
- Test: `tests/test_caom_shapes.py`

**Interfaces:**
- Consumes: `require_caom2` from Task 2.
- Produces:
  - `region_to_caom_shape(region: Region | None) -> Any | None` returning a `caom2.shape.Polygon`, `Circle`, `Box`, or `None`.
  - `arcsec_to_degrees(value: float | None) -> float | None`.
  - `metres(quantity: astropy.units.Quantity) -> float`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_caom_shapes.py`, and add these imports at the top of the file:

```python
import astropy.units as u
import lsst.sphgeom
```

```python
class ConversionTestCase(unittest.TestCase):
    """Tests of region and unit conversion."""

    def test_convex_polygon(self):
        """A ConvexPolygon becomes a caom2 Polygon with matching vertices."""
        vertices = [
            lsst.sphgeom.UnitVector3d(lsst.sphgeom.LonLat.fromDegrees(lon, lat))
            for lon, lat in ((10.0, 20.0), (11.0, 20.0), (11.0, 21.0), (10.0, 21.0))
        ]
        region = lsst.sphgeom.ConvexPolygon(vertices)

        shape = _caom_shapes.region_to_caom_shape(region)

        self.assertEqual(len(shape.points), 4)
        self.assertEqual(len(shape.samples.vertices), 5)
        expected = {(round(v.getLon().asDegrees(), 6), round(v.getLat().asDegrees(), 6))
                    for v in (lsst.sphgeom.LonLat(vec) for vec in region.getVertices())}
        actual = {(round(p.cval1, 6), round(p.cval2, 6)) for p in shape.points}
        self.assertEqual(actual, expected)

    def test_circle(self):
        """A sphgeom Circle becomes a caom2 Circle in degrees."""
        center = lsst.sphgeom.UnitVector3d(lsst.sphgeom.LonLat.fromDegrees(30.0, -10.0))
        region = lsst.sphgeom.Circle(center, lsst.sphgeom.Angle.fromDegrees(2.5))

        shape = _caom_shapes.region_to_caom_shape(region)

        self.assertAlmostEqual(shape.center.cval1, 30.0, places=6)
        self.assertAlmostEqual(shape.center.cval2, -10.0, places=6)
        self.assertAlmostEqual(shape.radius, 2.5, places=6)

    def test_none_region(self):
        """A missing region converts to None rather than raising."""
        self.assertIsNone(_caom_shapes.region_to_caom_shape(None))

    def test_unsupported_region_warns(self):
        """An unsupported region logs a warning and yields None."""
        a = lsst.sphgeom.Circle(
            lsst.sphgeom.UnitVector3d(lsst.sphgeom.LonLat.fromDegrees(0.0, 0.0)),
            lsst.sphgeom.Angle.fromDegrees(1.0),
        )
        b = lsst.sphgeom.Circle(
            lsst.sphgeom.UnitVector3d(lsst.sphgeom.LonLat.fromDegrees(90.0, 0.0)),
            lsst.sphgeom.Angle.fromDegrees(1.0),
        )
        region = lsst.sphgeom.UnionRegion(a, b)

        with self.assertLogs("lsst.dax.obscore._caom_shapes", level="WARNING"):
            self.assertIsNone(_caom_shapes.region_to_caom_shape(region))

    def test_arcsec_to_degrees(self):
        """Pixel scale converts from arcsec to degrees via astropy."""
        self.assertAlmostEqual(
            _caom_shapes.arcsec_to_degrees(0.2), (0.2 * u.arcsec).to_value(u.deg), places=12
        )
        self.assertIsNone(_caom_shapes.arcsec_to_degrees(None))
```

- [ ] **Step 2: Run to verify they fail**

Run: `pytest tests/test_caom_shapes.py -q`

Expected: FAIL with `AttributeError: module ... has no attribute 'region_to_caom_shape'`.

- [ ] **Step 3: Implement the helpers**

In `_caom_shapes.py`, extend `__all__` to `["arcsec_to_degrees", "metres", "region_to_caom_shape", "require_caom2"]` (ruff rule RUF022 requires it sorted), add the imports, and append the implementations:

```python
import astropy.units as u
from lsst.sphgeom import Box, Circle, ConvexPolygon, LonLat, Region
from lsst.utils.logging import getLogger

_LOG = getLogger(__name__)
```

```python
def arcsec_to_degrees(value: float | None) -> float | None:
    """Convert an angle in arcseconds to degrees.

    Parameters
    ----------
    value : `float` or `None`
        Angle in arcseconds. `None` passes through unchanged.

    Returns
    -------
    degrees : `float` or `None`
        The angle in degrees.
    """
    if value is None:
        return None
    return float((value * u.arcsec).to_value(u.deg))


def metres(quantity: u.Quantity) -> float:
    """Convert a length quantity to a plain number of metres.

    Parameters
    ----------
    quantity : `astropy.units.Quantity`
        A quantity convertible to metres.

    Returns
    -------
    value : `float`
        The value in metres.
    """
    return float(quantity.to_value(u.m))


def region_to_caom_shape(region: Region | None) -> Any | None:
    """Convert a sphgeom region to the equivalent CAOM shape.

    Parameters
    ----------
    region : `~lsst.sphgeom.Region` or `None`
        Region to convert. `None` passes through unchanged.

    Returns
    -------
    shape : `caom2.shape.Polygon`, `caom2.shape.Circle`, \
            `caom2.shape.Box` or `None`
        The equivalent CAOM shape, or `None` if the region is `None` or of
        an unsupported type.

    Notes
    -----
    An unsupported region type, such as a union region, issues a warning
    and returns `None` rather than raising, so that a single awkward
    dataset does not abort an entire export.
    """
    if region is None:
        return None
    require_caom2()

    if isinstance(region, ConvexPolygon):
        points = [
            caom2.shape.Point(coord.getLon().asDegrees(), coord.getLat().asDegrees())
            for coord in (LonLat(vector) for vector in region.getVertices())
        ]
        # A caom2 Polygon carries both a simple point list and a sampled
        # MultiPolygon; the latter is explicitly closed, so the first
        # vertex is repeated as the final LINE segment.
        vertices = [
            caom2.shape.Vertex(point.cval1, point.cval2, caom2.shape.SegmentType.LINE)
            for point in points
        ]
        vertices[0].type = caom2.shape.SegmentType.MOVE
        vertices.append(
            caom2.shape.Vertex(points[0].cval1, points[0].cval2, caom2.shape.SegmentType.CLOSE)
        )
        return caom2.shape.Polygon(points=points, samples=caom2.shape.MultiPolygon(vertices=vertices))

    if isinstance(region, Circle):
        center = LonLat(region.getCenter())
        return caom2.shape.Circle(
            center=caom2.shape.Point(center.getLon().asDegrees(), center.getLat().asDegrees()),
            radius=region.getOpeningAngle().asDegrees(),
        )

    if isinstance(region, Box):
        lon = region.getLon()
        lat = region.getLat()
        center = LonLat(region.getCenter())
        return caom2.shape.Box(
            center=caom2.shape.Point(center.getLon().asDegrees(), center.getLat().asDegrees()),
            width=(lon.getB() - lon.getA()).asDegrees(),
            height=(lat.getB() - lat.getA()).asDegrees(),
        )

    _LOG.warning("Cannot convert region of type %s to a CAOM shape; omitting position.", type(region))
    return None
```

Add `from typing import Any` to the imports.

- [ ] **Step 4: Run the tests to verify they pass**

Run: `pytest tests/test_caom_shapes.py -q`

Expected: 7 passed. If `test_convex_polygon` fails on the sample count, check that the closing vertex is appended exactly once.

- [ ] **Step 5: Confirm no bare conversion constants**

Run: `grep -nE '3600|180 */ *math\.pi|math\.pi */ *180' python/lsst/dax/obscore/_caom_shapes.py`

Expected: no output.

- [ ] **Step 6: Lint and commit**

Run: `ruff check python/lsst/dax/obscore/_caom_shapes.py tests/test_caom_shapes.py && ruff format --check python/lsst/dax/obscore/_caom_shapes.py && mypy python/lsst/dax/obscore/_caom_shapes.py`

```bash
git add python/lsst/dax/obscore/_caom_shapes.py tests/test_caom_shapes.py
git commit -m "Add sphgeom to CAOM shape and unit conversion helpers"
```

---

### Task 4: The `caom` configuration models

**Files:**
- Create: `python/lsst/dax/obscore/caom_config.py`
- Modify: `python/lsst/dax/obscore/config.py` (add the `caom` field to `ExporterConfig`)
- Test: `tests/test_caom_config.py`

**Interfaces:**
- Consumes: nothing. Must not import `caom2` or `_caom_shapes`.
- Produces:
  - `CaomDatasetTypeConfig` with fields `observation_id_fmt: str`, `product_id_fmt: str`, `algorithm: str | None`, `derived: bool`, `content_type: str | None`, `s_pixel_scale: float | None`, `observation_type_fmt: str | None`, `intent_fmt: str | None`, `auxiliary_datasets: dict[str, str]`.
  - `CaomProvenanceConfig` with `name: str`, `version: str | None`, `project: str | None`.
  - `CaomConfig` with `telescope_name: str | None`, `geo_location: tuple[float, float, float] | None`, `proposal_id: str | None`, `em_band: str | None`, `provenance: CaomProvenanceConfig | None`, `read_groups: list[str]`, `intent: str`, `observation_type: str`, `content_type: str | None`, `artifact_uri_fmt: str | None`, `dataset_types: dict[str, CaomDatasetTypeConfig]`, plus method `pixel_scale_for(dataset_type: str) -> float | None`.
  - `ExporterConfig.caom: CaomConfig | None`.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_caom_config.py` with the GPL header, then:

```python
import unittest

import pydantic

from lsst.dax.obscore import ExporterConfig
from lsst.dax.obscore.caom_config import CaomConfig
from lsst.dax.obscore.tests import DaxObsCoreTestMixin


class CaomConfigTestCase(unittest.TestCase, DaxObsCoreTestMixin):
    """Tests of the caom configuration block."""

    def test_absent_by_default(self):
        """Configurations without a caom block are still valid."""
        config = self.make_export_config()
        self.assertIsNone(config.caom)

    def test_unknown_dataset_type_rejected(self):
        """A caom dataset type must exist in the ObsCore dataset types."""
        config = self.make_export_config()
        with self.assertRaises(pydantic.ValidationError) as cm:
            ExporterConfig.model_validate(
                config.model_dump()
                | {
                    "caom": {
                        "dataset_types": {
                            "not_a_dataset_type": {
                                "observation_id_fmt": "{visit}",
                                "product_id_fmt": "{detector}",
                            }
                        }
                    }
                }
            )
        self.assertIn("not_a_dataset_type", str(cm.exception))

    def test_algorithm_requires_derived(self):
        """An algorithm on a non-derived dataset type is a config error."""
        with self.assertRaises(pydantic.ValidationError):
            CaomConfig.model_validate(
                {
                    "dataset_types": {
                        "deep_coadd": {
                            "observation_id_fmt": "{skymap}-{tract}",
                            "product_id_fmt": "{patch}-{band}",
                            "algorithm": "some.algorithm",
                        }
                    }
                }
            )

    def test_defaults(self):
        """Unset optional fields take the documented defaults."""
        config = CaomConfig.model_validate(
            {
                "dataset_types": {
                    "deep_coadd": {
                        "observation_id_fmt": "{skymap}-{tract}",
                        "product_id_fmt": "{patch}-{band}",
                    }
                }
            }
        )
        self.assertEqual(config.intent, "science")
        self.assertEqual(config.observation_type, "science")
        self.assertEqual(config.read_groups, [])
        self.assertIsNone(config.telescope_name)
        self.assertFalse(config.dataset_types["deep_coadd"].derived)

    def test_pixel_scale_for(self):
        """The placeholder pixel scale is read through one accessor."""
        config = CaomConfig.model_validate(
            {
                "dataset_types": {
                    "deep_coadd": {
                        "observation_id_fmt": "{skymap}-{tract}",
                        "product_id_fmt": "{patch}-{band}",
                        "s_pixel_scale": 0.2,
                    },
                    "visit_image": {
                        "observation_id_fmt": "{visit}",
                        "product_id_fmt": "{detector}",
                    },
                }
            }
        )
        self.assertEqual(config.pixel_scale_for("deep_coadd"), 0.2)
        self.assertIsNone(config.pixel_scale_for("visit_image"))


if __name__ == "__main__":
    unittest.main()
```

- [ ] **Step 2: Run to verify they fail**

Run: `pytest tests/test_caom_config.py -q`

Expected: FAIL with `ModuleNotFoundError: No module named 'lsst.dax.obscore.caom_config'`.

- [ ] **Step 3: Create `caom_config.py`**

Create the file with the GPL header, then:

```python
from __future__ import annotations

__all__ = ["CaomConfig", "CaomDatasetTypeConfig", "CaomProvenanceConfig"]

from pydantic import BaseModel, Field, model_validator


class CaomProvenanceConfig(BaseModel):
    """Provenance values that apply to every plane in the export."""

    name: str
    """Collection-specific common name of the process that produced the
    data, for example ``LSST Science Pipelines``. This names the software,
    not the dataset."""

    version: str | None = None
    """Version of that software, for example ``v29.0``. This is a fixed
    release-level value; the per-execution value is the Butler RUN
    collection, which is reported as ``Provenance.runID``."""

    project: str | None = None
    """Project the process belongs to, for example ``DRP``."""


class CaomDatasetTypeConfig(BaseModel):
    """CAOM configuration for a single dataset type."""

    observation_id_fmt: str
    """Template for the CAOM ``Observation.observationID``. Dataset types
    whose templates expand to the same value contribute planes to the same
    Observation."""

    product_id_fmt: str
    """Template for the CAOM ``Plane.productID``, unique within an
    Observation."""

    algorithm: str | None = None
    """Algorithm name for a derived Observation. Only valid when
    ``derived`` is true."""

    derived: bool = False
    """If true, emit a ``DerivedObservation`` rather than a
    ``SimpleObservation``."""

    content_type: str | None = None
    """MIME type of the files, overriding the top-level default."""

    s_pixel_scale: float | None = None
    """Pixel scale in arcsec, used for ``Position.sampleSize``.

    This is a placeholder. The value belongs beside ``s_xel`` in the
    ObsCore ``DatasetTypeConfig``, and moves there once daf_butler
    supports it. Read it through `CaomConfig.pixel_scale_for` so that the
    migration touches one function.
    """

    observation_type_fmt: str | None = None
    """Template for ``Observation.type``, overriding the top-level
    literal. Use this for exposure-based dataset types, where the value
    is carried on the exposure dimension record."""

    intent_fmt: str | None = None
    """Template for ``Observation.intent``, overriding the top-level
    literal. Must expand to ``science`` or ``calibration``."""

    auxiliary_datasets: dict[str, str] = Field(default_factory=dict)
    """Dataset types contributing auxiliary artifacts to each plane.

    Keys are dataset type names. Values must be member values of
    `caom2.ProductType`, which does not include ``science``; use
    ``auxiliary``, ``weight``, ``preview`` and so on. The plane's primary
    artifact always uses ``this``.
    """

    @model_validator(mode="after")
    def _check_algorithm(self) -> CaomDatasetTypeConfig:
        """Reject an algorithm on a non-derived dataset type.

        Returns
        -------
        config : `CaomDatasetTypeConfig`
            This configuration, unchanged.
        """
        if self.algorithm is not None and not self.derived:
            raise ValueError(
                "'algorithm' is only meaningful when 'derived' is true; "
                "caom2 forces a SimpleObservation to use the algorithm name 'exposure'."
            )
        if self.derived and self.algorithm is None:
            raise ValueError("'derived' dataset types must specify an 'algorithm'.")
        return self


class CaomConfig(BaseModel):
    """CAOM extensions to an ObsCore exporter configuration.

    Notes
    -----
    Values already present in the ObsCore configuration are not repeated
    here. ``obs_collection`` supplies the CAOM collection,
    ``facility_name`` the telescope name, and the ``origin`` block the
    proposal title, producer, provenance reference and metadata release
    date.
    """

    telescope_name: str | None = None
    """CAOM ``Telescope.name``. Defaults to the ObsCore facility name."""

    geo_location: tuple[float, float, float] | None = None
    """Geocentric telescope position in metres, overriding the lookup by
    facility name. Needed only for facilities astropy does not know."""

    proposal_id: str | None = None
    """CAOM ``Proposal.id``."""

    em_band: str | None = None
    """CAOM ``EnergyBand`` name, for example ``OPTICAL``."""

    provenance: CaomProvenanceConfig | None = None
    """Provenance values shared by every plane."""

    read_groups: list[str] = Field(default_factory=list)
    """Group URIs controlling access to the data."""

    intent: str = "science"
    """Default ``Observation.intent``."""

    observation_type: str = "science"
    """Default ``Observation.type``.

    This matches the vocabulary of the ``exposure`` dimension's
    ``observation_type`` field, so that literal and templated values agree.
    """

    content_type: str | None = None
    """Default MIME type of the exported files."""

    artifact_uri_fmt: str | None = None
    """Template for ``Artifact.uri``. ``{butler_uri}`` is available when a
    datastore is present. If unset, the resolved Butler URI is used."""

    dataset_types: dict[str, CaomDatasetTypeConfig]
    """Per-dataset-type configuration. Only these dataset types are
    exported."""

    def pixel_scale_for(self, dataset_type: str) -> float | None:
        """Return the pixel scale in arcsec for a dataset type.

        Parameters
        ----------
        dataset_type : `str`
            Name of the dataset type.

        Returns
        -------
        pixel_scale : `float` or `None`
            Pixel scale in arcsec, or `None` if not configured.

        Notes
        -----
        This is the single point of access for the placeholder
        ``s_pixel_scale`` key. When the value moves to the ObsCore
        ``DatasetTypeConfig``, only this method changes.
        """
        config = self.dataset_types.get(dataset_type)
        return config.s_pixel_scale if config is not None else None
```

- [ ] **Step 4: Wire it into `ExporterConfig`**

In `python/lsst/dax/obscore/config.py`, add `from .caom_config import CaomConfig` to the imports, and add to `ExporterConfig` immediately after the `origin` field:

```python
    caom: CaomConfig | None = None
    """Optional CAOM extensions, used only by the CAOM exporter."""

    @model_validator(mode="after")
    def _check_caom_dataset_types(self) -> ExporterConfig:
        """Check CAOM dataset types are known to the ObsCore config.

        Returns
        -------
        config : `ExporterConfig`
            This configuration, unchanged.
        """
        if self.caom is not None:
            unknown = set(self.caom.dataset_types) - set(self.dataset_types)
            if unknown:
                raise ValueError(
                    f"CAOM dataset types {sorted(unknown)} are not defined in 'dataset_types'."
                )
        return self
```

Add `model_validator` to the existing pydantic import line: `from pydantic import BaseModel, ConfigDict, Field, model_validator`.

- [ ] **Step 5: Run the tests to verify they pass**

Run: `pytest tests/test_caom_config.py -q`

Expected: 5 passed.

- [ ] **Step 6: Check nothing else broke and that `caom2` is not imported**

Run: `pytest tests/ -q`

Expected: all pass.

Run: `python -c "import lsst.dax.obscore.caom_config, sys; assert 'caom2' not in sys.modules; print('clean')"`

Expected: `clean`.

- [ ] **Step 7: Lint and commit**

Run: `ruff check python/lsst/dax/obscore/caom_config.py python/lsst/dax/obscore/config.py tests/test_caom_config.py && mypy python/lsst/dax/obscore/caom_config.py python/lsst/dax/obscore/config.py`

```bash
git add python/lsst/dax/obscore/caom_config.py python/lsst/dax/obscore/config.py tests/test_caom_config.py
git commit -m "Add caom configuration block to ExporterConfig"
```

---

### Task 5: Build Observations for a single dataset type

The core of the exporter. Grouping across dataset types comes in Task 6; this task handles one dataset type at a time.

**Files:**
- Create: `python/lsst/dax/obscore/caom_exporter.py`
- Modify: `python/lsst/dax/obscore/tests.py` (add `make_caom_config`)
- Test: `tests/test_caom_exporter.py`

**Interfaces:**
- Consumes: `ObscoreExporter._iter_record_refs` and `_QueryState` (Task 1); `region_to_caom_shape`, `arcsec_to_degrees`, `metres`, `require_caom2` (Tasks 2-3); `CaomConfig` (Task 4).
- Produces: `CaomExporter(butler: Butler, config: ExporterConfig)` with `iter_observations() -> Iterator[Any]` yielding `caom2.Observation` instances.

- [ ] **Step 1: Add the test configuration helper**

In `python/lsst/dax/obscore/tests.py`, add this method to `DaxObsCoreTestMixin`:

```python
    def make_caom_config(self) -> ExporterConfig:
        """Return an exporter configuration with a CAOM block.

        Returns
        -------
        config : `ExporterConfig`
            Configuration covering the mock dataset types, extended with
            CAOM settings.
        """
        config = self.make_export_config()
        config.caom = CaomConfig.model_validate(
            {
                "telescope_name": "Subaru Telescope",
                "proposal_id": "TEST-PROPOSAL",
                "em_band": "OPTICAL",
                "provenance": {"name": "LSST Science Pipelines", "version": "v29.0"},
                "read_groups": ["ivo://example.org/gms?TEST"],
                "content_type": "application/fits",
                "artifact_uri_fmt": "cadc:TEST/{dataset_type}/{id}",
                "dataset_types": {
                    "_mock_calexp": {
                        "observation_id_fmt": "{records[visit].name}",
                        "product_id_fmt": "calexp-{detector}",
                        "s_pixel_scale": 0.17,
                    },
                    "_mock_deepCoadd": {
                        "observation_id_fmt": "{skymap}-{tract}",
                        "product_id_fmt": "deepCoadd-{patch}-{band}",
                        "derived": True,
                        "algorithm": "mock.coadd",
                        "s_pixel_scale": 0.17,
                    },
                },
            }
        )
        return config
```

Add `from .caom_config import CaomConfig` to the imports in that file.

- [ ] **Step 2: Write the failing tests**

Create `tests/test_caom_exporter.py` with the GPL header, then:

```python
import os
import unittest

from lsst.daf.butler.tests.utils import makeTestTempDir, removeTestTempDir
from lsst.dax.obscore.tests import DaxObsCoreTestMixin

try:
    import caom2
except ImportError:
    caom2 = None

TESTDIR = os.path.abspath(os.path.dirname(__file__))


@unittest.skipUnless(caom2 is not None, "caom2 is not installed")
class CaomExporterTestCase(unittest.TestCase, DaxObsCoreTestMixin):
    """Tests of CaomExporter."""

    def setUp(self):
        self.root = makeTestTempDir(TESTDIR)

    def tearDown(self):
        removeTestTempDir(self.root)

    def make_populated_butler(self):
        """Return a butler with the test data imported."""
        butler = self.make_butler()
        self.enterContext(butler)
        butler.import_(
            filename=os.path.join(TESTDIR, "data", "hsc_gen3.yaml"), without_datastore=True
        )
        return butler

    def test_calexp_observations(self):
        """Visit-based records group into one Observation per visit."""
        from lsst.dax.obscore.caom_exporter import CaomExporter

        butler = self.make_populated_butler()
        config = self.make_caom_config()
        config.select_dataset_types(["_mock_calexp"])
        config.caom.dataset_types = {"_mock_calexp": config.caom.dataset_types["_mock_calexp"]}

        observations = list(CaomExporter(butler, config).iter_observations())

        self.assertGreater(len(observations), 0)
        for observation in observations:
            self.assertIsInstance(observation, caom2.SimpleObservation)
            self.assertEqual(observation.collection, "obs-collection")
            self.assertEqual(observation.telescope.name, "Subaru Telescope")
            self.assertEqual(observation.proposal.id, "TEST-PROPOSAL")
            self.assertEqual(observation.type, "science")
            self.assertEqual(observation.intent, caom2.ObservationIntentType.SCIENCE)
            self.assertGreater(len(observation.planes), 0)
            for plane in observation.planes.values():
                self.assertEqual(plane.calibration_level, caom2.CalibrationLevel.CALIBRATED)
                self.assertEqual(plane.data_product_type, caom2.DataProductType.IMAGE)
                self.assertEqual(plane.provenance.name, "LSST Science Pipelines")
                self.assertEqual(plane.provenance.version, "v29.0")
                self.assertIsNotNone(plane.provenance.run_id)
                self.assertIsNotNone(plane.position)
                self.assertIsInstance(plane.position.bounds, caom2.shape.Polygon)
                self.assertIsNotNone(plane.energy)
                self.assertIsNotNone(plane.time)
                self.assertGreaterEqual(len(plane.artifacts), 1)
                for artifact in plane.artifacts.values():
                    self.assertTrue(artifact.uri.startswith("cadc:TEST/_mock_calexp/"))
                    self.assertEqual(artifact.content_type, "application/fits")
                    self.assertEqual(artifact.release_type, caom2.ReleaseType.DATA)
                    self.assertEqual(artifact.product_type, caom2.ProductType.THIS)

    def test_coadd_is_derived(self):
        """Coadd records produce DerivedObservations keyed on tract."""
        from lsst.dax.obscore.caom_exporter import CaomExporter

        butler = self.make_populated_butler()
        config = self.make_caom_config()
        config.select_dataset_types(["_mock_deepCoadd"])
        config.caom.dataset_types = {"_mock_deepCoadd": config.caom.dataset_types["_mock_deepCoadd"]}

        observations = list(CaomExporter(butler, config).iter_observations())

        self.assertGreater(len(observations), 0)
        for observation in observations:
            self.assertIsInstance(observation, caom2.DerivedObservation)
            self.assertEqual(observation.algorithm.name, "mock.coadd")
            # observation_id is skymap-tract, so it has exactly one hyphen
            # more than the skymap name contains.
            self.assertTrue(observation.observation_id.startswith("skymap-"))

    def test_pixel_scale_becomes_sample_size(self):
        """The configured pixel scale reaches Position.sampleSize in deg."""
        import astropy.units as u

        from lsst.dax.obscore.caom_exporter import CaomExporter

        butler = self.make_populated_butler()
        config = self.make_caom_config()
        config.select_dataset_types(["_mock_calexp"])
        config.caom.dataset_types = {"_mock_calexp": config.caom.dataset_types["_mock_calexp"]}

        observations = list(CaomExporter(butler, config).iter_observations())
        plane = next(iter(observations[0].planes.values()))

        self.assertAlmostEqual(
            plane.position.sample_size, (0.17 * u.arcsec).to_value(u.deg), places=12
        )


if __name__ == "__main__":
    unittest.main()
```

- [ ] **Step 3: Run to verify they fail**

Run: `pytest tests/test_caom_exporter.py -q`

Expected: FAIL with `ModuleNotFoundError: No module named 'lsst.dax.obscore.caom_exporter'`.

- [ ] **Step 4: Implement `CaomExporter`**

Create `python/lsst/dax/obscore/caom_exporter.py` with the GPL header, then:

```python
from __future__ import annotations

__all__ = ["CaomExporter"]

import datetime
from collections.abc import Iterator
from typing import Any

import astropy.units as u
from astropy.coordinates import EarthLocation

from lsst.daf.butler import Butler, DatasetRef
from lsst.sphgeom import Region
from lsst.utils.logging import getLogger

from ._caom_shapes import arcsec_to_degrees, caom2, metres, region_to_caom_shape, require_caom2
from .config import ExporterConfig
from .obscore_exporter import ObscoreExporter, _QueryState

_LOG = getLogger(__name__)


class CaomExporter:
    """Export Butler datasets as CAOM Observations.

    Parameters
    ----------
    butler : `lsst.daf.butler.Butler`
        Data butler.
    config : `lsst.dax.obscore.ExporterConfig`
        Exporter configuration. Must have a ``caom`` block.

    Raises
    ------
    ValueError
        Raised if the configuration has no ``caom`` block.
    ImportError
        Raised if the optional ``caom2`` dependency is not installed.
    """

    def __init__(self, butler: Butler, config: ExporterConfig):
        require_caom2()
        if config.caom is None:
            raise ValueError("Configuration has no 'caom' section; cannot export CAOM.")

        self.butler = butler
        self.config = config
        self.caom_config = config.caom

        # Only export dataset types with CAOM configuration.
        config.select_dataset_types(self.caom_config.dataset_types)
        self._obscore = ObscoreExporter(butler, config)

    def iter_observations(self) -> Iterator[Any]:
        """Generate CAOM Observations for the configured dataset types.

        Yields
        ------
        observation : `caom2.Observation`
            One Observation per distinct expanded ``observation_id_fmt``.
        """
        observations: dict[str, Any] = {}
        state = _QueryState()
        for ref, region, record in self._obscore._iter_record_refs(state):
            self._add_record(observations, ref, region, record)
        yield from observations.values()

    def _format_keywords(self, ref: DatasetRef, record: dict[str, Any]) -> dict[str, Any]:
        """Build the namespace used to expand configuration templates.

        Parameters
        ----------
        ref : `~lsst.daf.butler.DatasetRef`
            Reference to the dataset.
        record : `dict` [ `str`, `~typing.Any` ]
            The completed ObsCore record.

        Returns
        -------
        keywords : `dict` [ `str`, `~typing.Any` ]
            Values available to every template, comprising the data ID
            mapping, its dimension records, the dataset identity, and every
            ObsCore column.
        """
        keywords: dict[str, Any] = dict(records=ref.dataId.records)
        keywords.update(ref.dataId.mapping)
        keywords.update(id=ref.id, run=ref.run, dataset_type=ref.datasetType.name)
        keywords.update(record)
        return keywords

    def _add_record(
        self,
        observations: dict[str, Any],
        ref: DatasetRef,
        region: Region | None,
        record: dict[str, Any],
    ) -> None:
        """Fold one ObsCore record into the accumulating Observations.

        Parameters
        ----------
        observations : `dict` [ `str`, `caom2.Observation` ]
            Observations accumulated so far, keyed by observation ID.
            Updated in place.
        ref : `~lsst.daf.butler.DatasetRef`
            Reference to the dataset.
        region : `~lsst.sphgeom.Region` or `None`
            Spatial region for the dataset.
        record : `dict` [ `str`, `~typing.Any` ]
            The completed ObsCore record.
        """
        dataset_type = ref.datasetType.name
        dataset_config = self.caom_config.dataset_types[dataset_type]
        keywords = self._format_keywords(ref, record)

        observation_id = dataset_config.observation_id_fmt.format(**keywords)
        product_id = dataset_config.product_id_fmt.format(**keywords)

        observation = observations.get(observation_id)
        if observation is None:
            observation = self._make_observation(observation_id, dataset_config, keywords, record)
            observations[observation_id] = observation

        plane = observation.planes.get(product_id)
        if plane is None:
            plane = self._make_plane(product_id, dataset_config, ref, region, record)
            observation.planes[product_id] = plane

        artifact = self._make_artifact(dataset_config, keywords, product_type_name="this")
        if artifact is not None:
            plane.artifacts[artifact.uri] = artifact

    def _make_observation(
        self,
        observation_id: str,
        dataset_config: Any,
        keywords: dict[str, Any],
        record: dict[str, Any],
    ) -> Any:
        """Create a CAOM Observation.

        Parameters
        ----------
        observation_id : `str`
            Expanded observation identifier.
        dataset_config : `CaomDatasetTypeConfig`
            CAOM configuration for the dataset type of the first record
            encountered for this observation.
        keywords : `dict` [ `str`, `~typing.Any` ]
            Template namespace for this record.
        record : `dict` [ `str`, `~typing.Any` ]
            The completed ObsCore record.

        Returns
        -------
        observation : `caom2.Observation`
            A new Observation with no planes.
        """
        config = self.caom_config
        observation_type = (
            dataset_config.observation_type_fmt.format(**keywords)
            if dataset_config.observation_type_fmt
            else config.observation_type
        )
        intent_name = (
            dataset_config.intent_fmt.format(**keywords) if dataset_config.intent_fmt else config.intent
        )
        intent = caom2.ObservationIntentType(intent_name.lower())

        target = None
        if record.get("target_name"):
            target = caom2.Target(name=record["target_name"])

        instrument = None
        if record.get("instrument_name"):
            instrument = caom2.Instrument(record["instrument_name"])

        proposal = None
        if config.proposal_id is not None:
            origin = self.config.origin
            proposal = caom2.Proposal(
                id=config.proposal_id,
                project=self.config.obs_collection,
                title=origin.title if origin is not None else None,
            )

        meta_release = None
        if self.config.origin is not None:
            meta_release = datetime.datetime.combine(
                self.config.origin.publication_date, datetime.time(), tzinfo=datetime.UTC
            )

        kwargs: dict[str, Any] = dict(
            collection=self.config.obs_collection,
            observation_id=observation_id,
            intent=intent,
            type=observation_type,
            proposal=proposal,
            telescope=self._make_telescope(record),
            instrument=instrument,
            target=target,
            meta_release=meta_release,
        )
        if dataset_config.derived:
            return caom2.DerivedObservation(algorithm=caom2.Algorithm(dataset_config.algorithm), **kwargs)
        return caom2.SimpleObservation(**kwargs)

    def _make_telescope(self, record: dict[str, Any]) -> Any:
        """Create a CAOM Telescope with a geocentric position.

        Parameters
        ----------
        record : `dict` [ `str`, `~typing.Any` ]
            The completed ObsCore record, used for the facility name.

        Returns
        -------
        telescope : `caom2.Telescope`
            The telescope, with geocentric coordinates when they can be
            resolved.
        """
        facility = record.get("facility_name") or self.config.facility_name
        name = self.caom_config.telescope_name or facility

        location = self.caom_config.geo_location
        if location is None:
            try:
                site = EarthLocation.of_site(facility)
            except Exception:
                _LOG.warning(
                    "Could not resolve a geocentric position for facility %r; "
                    "set caom.geo_location to supply one.",
                    facility,
                )
                return caom2.Telescope(name=name)
            location = tuple(metres(value) for value in site.geocentric)

        return caom2.Telescope(
            name=name,
            geo_location_x=location[0],
            geo_location_y=location[1],
            geo_location_z=location[2],
        )

    def _make_plane(
        self,
        product_id: str,
        dataset_config: Any,
        ref: DatasetRef,
        region: Region | None,
        record: dict[str, Any],
    ) -> Any:
        """Create a CAOM Plane with its position, energy and time.

        Parameters
        ----------
        product_id : `str`
            Expanded product identifier.
        dataset_config : `CaomDatasetTypeConfig`
            CAOM configuration for the dataset type.
        ref : `~lsst.daf.butler.DatasetRef`
            Reference to the dataset.
        region : `~lsst.sphgeom.Region` or `None`
            Spatial region for the dataset.
        record : `dict` [ `str`, `~typing.Any` ]
            The completed ObsCore record.

        Returns
        -------
        plane : `caom2.Plane`
            A new plane with no artifacts.
        """
        config = self.caom_config
        origin = self.config.origin

        provenance = None
        if config.provenance is not None:
            provenance = caom2.Provenance(
                name=config.provenance.name,
                version=config.provenance.version,
                project=config.provenance.project,
                producer=origin.publisher if origin is not None else None,
                run_id=ref.run,
                reference=origin.reference_url if origin is not None else None,
            )

        observable = None
        if record.get("o_ucd"):
            observable = caom2.Observable(record["o_ucd"])

        read_groups = caom2.caom_util.URISet()
        for uri in config.read_groups:
            read_groups.add(uri)

        plane = caom2.Plane(
            product_id=product_id,
            data_product_type=caom2.DataProductType(record["dataproduct_type"]),
            calibration_level=caom2.CalibrationLevel(record["calib_level"]),
            provenance=provenance,
            observable=observable,
            data_read_groups=read_groups,
        )

        bounds = region_to_caom_shape(region)
        dimension = None
        if record.get("s_xel1") is not None and record.get("s_xel2") is not None:
            dimension = caom2.wcs.Dimension2D(record["s_xel1"], record["s_xel2"])
        sample_size = arcsec_to_degrees(config.pixel_scale_for(ref.datasetType.name))
        if bounds is not None or dimension is not None or sample_size is not None:
            plane.position = caom2.Position(
                bounds=bounds,
                dimension=dimension,
                resolution=record.get("s_resolution"),
                sample_size=sample_size,
            )

        em_min = record.get("em_min")
        em_max = record.get("em_max")
        if em_min is not None and em_max is not None:
            # ObsCore em_min and em_max are already metres, as CAOM
            # expects, but the conversion is stated rather than assumed.
            lower = metres(em_min * u.m)
            upper = metres(em_max * u.m)
            resolving_power = ((upper + lower) / 2.0) / (upper - lower) if upper > lower else None
            plane.energy = caom2.Energy(
                bounds=caom2.shape.Interval(lower, upper),
                bandpass_name=record.get("em_filter_name"),
                resolving_power=resolving_power,
                em_band=caom2.EnergyBand[config.em_band] if config.em_band else None,
            )

        t_min = record.get("t_min")
        t_max = record.get("t_max")
        if t_min is not None and t_max is not None:
            plane.time = caom2.Time(
                bounds=caom2.shape.Interval(t_min, t_max),
                exposure=record.get("t_exptime"),
            )

        return plane

    def _make_artifact(
        self, dataset_config: Any, keywords: dict[str, Any], product_type_name: str
    ) -> Any | None:
        """Create a CAOM Artifact for one dataset.

        Parameters
        ----------
        dataset_config : `CaomDatasetTypeConfig`
            CAOM configuration for the dataset type.
        keywords : `dict` [ `str`, `~typing.Any` ]
            Template namespace for this record.
        product_type_name : `str`
            CAOM product type name. Must be a member value of
            `caom2.ProductType`, which has no ``science`` term: the primary
            artifact of a plane uses ``this`` and extras use ``auxiliary``.

        Returns
        -------
        artifact : `caom2.Artifact` or `None`
            The artifact, or `None` if no URI could be determined.
        """
        uri_fmt = self.caom_config.artifact_uri_fmt
        if uri_fmt is None:
            _LOG.warning(
                "No caom.artifact_uri_fmt configured and no Butler URI available for %s; "
                "skipping artifact.",
                keywords["id"],
            )
            return None
        uri = uri_fmt.format(**keywords)

        return caom2.Artifact(
            uri=uri,
            product_type=caom2.ProductType(product_type_name),
            release_type=caom2.ReleaseType.DATA,
            content_type=dataset_config.content_type or self.caom_config.content_type,
        )
```

- [ ] **Step 5: Run the tests to verify they pass**

Run: `pytest tests/test_caom_exporter.py -q`

Expected: 3 passed.

If `caom2.wcs.Dimension2D` is not importable at that path, run `python -c "import caom2; print(caom2.Dimension2D)"` and use whichever path resolves.

- [ ] **Step 6: Confirm nothing else broke**

Run: `pytest tests/ -q`

Expected: all pass.

- [ ] **Step 7: Lint and commit**

Run: `ruff check python/lsst/dax/obscore/caom_exporter.py python/lsst/dax/obscore/tests.py tests/test_caom_exporter.py && ruff format --check python/lsst/dax/obscore/caom_exporter.py && mypy python/lsst/dax/obscore/caom_exporter.py`

```bash
git add python/lsst/dax/obscore/caom_exporter.py python/lsst/dax/obscore/tests.py tests/test_caom_exporter.py
git commit -m "Add CaomExporter building Observations from ObsCore records"
```

---

### Task 6: Cross-dataset-type merging and conflict rules

Task 5 already merges by observation ID, because `iter_observations` keys a single dictionary. This task adds the conflict detection the spec requires and proves the merge works across dataset types.

**Files:**
- Modify: `python/lsst/dax/obscore/caom_exporter.py` (`_add_record`)
- Test: `tests/test_caom_exporter.py`

**Interfaces:**
- Consumes: `CaomExporter._add_record` from Task 5.
- Produces: no new public API. `_add_record` raises `ValueError` on a duplicate `(observation_id, product_id)` from a different dataset type, and warns on conflicting observation-level attributes.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_caom_exporter.py`, inside `CaomExporterTestCase`:

```python
    def test_merges_across_dataset_types(self):
        """Two dataset types sharing an observation ID share an Observation."""
        from lsst.dax.obscore.caom_exporter import CaomExporter

        butler = self.make_populated_butler()
        config = self.make_caom_config()
        # Give the coadd the same observation ID template as the calexp so
        # that the two dataset types collide deliberately.
        config.caom.dataset_types["_mock_deepCoadd"].observation_id_fmt = "shared-observation"
        config.caom.dataset_types["_mock_calexp"].observation_id_fmt = "shared-observation"

        observations = list(CaomExporter(butler, config).iter_observations())

        self.assertEqual(len(observations), 1)
        product_ids = set(observations[0].planes)
        self.assertTrue(any(p.startswith("calexp-") for p in product_ids))
        self.assertTrue(any(p.startswith("deepCoadd-") for p in product_ids))
        levels = {plane.calibration_level for plane in observations[0].planes.values()}
        self.assertEqual(
            levels, {caom2.CalibrationLevel.CALIBRATED, caom2.CalibrationLevel.PRODUCT}
        )

    def test_duplicate_product_id_is_an_error(self):
        """Colliding product IDs from different dataset types abort."""
        from lsst.dax.obscore.caom_exporter import CaomExporter

        butler = self.make_populated_butler()
        config = self.make_caom_config()
        config.caom.dataset_types["_mock_deepCoadd"].observation_id_fmt = "shared-observation"
        config.caom.dataset_types["_mock_calexp"].observation_id_fmt = "shared-observation"
        config.caom.dataset_types["_mock_deepCoadd"].product_id_fmt = "collide"
        config.caom.dataset_types["_mock_calexp"].product_id_fmt = "collide"

        with self.assertRaises(ValueError) as cm:
            list(CaomExporter(butler, config).iter_observations())
        self.assertIn("collide", str(cm.exception))
```

- [ ] **Step 2: Run to verify they fail**

Run: `pytest tests/test_caom_exporter.py -q -k "merges or duplicate"`

Expected: `test_merges_across_dataset_types` may pass already; `test_duplicate_product_id_is_an_error` FAILS because no `ValueError` is raised.

- [ ] **Step 3: Add conflict detection**

In `caom_exporter.py`, add an instance attribute in `__init__` after `self._obscore = ...`:

```python
        # Records which dataset type created each plane, so that a
        # product ID collision between dataset types can be reported.
        self._plane_owners: dict[tuple[str, str], str] = {}
```

Replace the plane lookup block in `_add_record` with:

```python
        owner_key = (observation_id, product_id)
        owner = self._plane_owners.get(owner_key)
        if owner is None:
            self._plane_owners[owner_key] = dataset_type
        elif owner != dataset_type:
            raise ValueError(
                f"Product ID {product_id!r} in observation {observation_id!r} is produced by both "
                f"{owner!r} and {dataset_type!r}. Dataset types sharing an observation must use "
                "distinct 'product_id_fmt' templates."
            )

        plane = observation.planes.get(product_id)
        if plane is None:
            plane = self._make_plane(product_id, dataset_config, ref, region, record)
            observation.planes[product_id] = plane
```

Immediately after the `observations[observation_id] = observation` branch, add the observation-level conflict warning:

```python
        else:
            self._check_observation_conflict(observation, dataset_config, keywords)
```

and add the method:

```python
    def _check_observation_conflict(
        self, observation: Any, dataset_config: Any, keywords: dict[str, Any]
    ) -> None:
        """Warn if a later dataset type disagrees on observation attributes.

        Parameters
        ----------
        observation : `caom2.Observation`
            The Observation created by an earlier dataset type.
        dataset_config : `CaomDatasetTypeConfig`
            CAOM configuration for the dataset type of the current record.
        keywords : `dict` [ `str`, `~typing.Any` ]
            Template namespace for the current record.

        Notes
        -----
        The first dataset type to reach an observation ID sets the
        observation-level attributes. A later disagreement is logged and
        ignored, because the prototype has no basis for choosing between
        them.
        """
        expected_derived = isinstance(observation, caom2.DerivedObservation)
        if expected_derived != dataset_config.derived:
            _LOG.warning(
                "Observation %s was created as %s but dataset type %s configures derived=%s; "
                "keeping the first.",
                observation.observation_id,
                type(observation).__name__,
                keywords["dataset_type"],
                dataset_config.derived,
            )

        if dataset_config.observation_type_fmt:
            observation_type = dataset_config.observation_type_fmt.format(**keywords)
            if observation_type != observation.type:
                _LOG.warning(
                    "Observation %s has type %r but dataset type %s gives %r; keeping the first.",
                    observation.observation_id,
                    observation.type,
                    keywords["dataset_type"],
                    observation_type,
                )
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `pytest tests/test_caom_exporter.py -q`

Expected: 5 passed.

- [ ] **Step 5: Lint and commit**

Run: `ruff check python/lsst/dax/obscore/caom_exporter.py tests/test_caom_exporter.py && mypy python/lsst/dax/obscore/caom_exporter.py`

```bash
git add python/lsst/dax/obscore/caom_exporter.py tests/test_caom_exporter.py
git commit -m "Detect product ID and observation attribute conflicts when merging"
```

---

### Task 7: Optional Butler URI resolution

Adds `{butler_uri}` to the template namespace when a datastore is available. The helper lives on `ObscoreExporter` so both exports can use it, per the spec.

**Files:**
- Modify: `python/lsst/dax/obscore/obscore_exporter.py` (add `_resolve_uris`)
- Modify: `python/lsst/dax/obscore/caom_exporter.py` (`iter_observations`, `_format_keywords`)
- Test: `tests/test_caom_exporter.py`

**Interfaces:**
- Consumes: `ObscoreExporter` from Task 1.
- Produces: `ObscoreExporter._resolve_uris(refs: list[DatasetRef]) -> dict[DatasetId, str]`, returning an empty mapping when no datastore is available.

- [ ] **Step 1: Write the failing test**

Append to `CaomExporterTestCase` in `tests/test_caom_exporter.py`:

```python
    def test_uri_template_without_datastore(self):
        """A datastore-less butler still exports, using the URI template."""
        from lsst.dax.obscore.caom_exporter import CaomExporter

        butler = self.make_populated_butler()
        config = self.make_caom_config()
        config.select_dataset_types(["_mock_calexp"])
        config.caom.dataset_types = {"_mock_calexp": config.caom.dataset_types["_mock_calexp"]}

        exporter = CaomExporter(butler, config)
        with self.assertLogs("lsst.dax.obscore.obscore_exporter", level="WARNING"):
            observations = list(exporter.iter_observations())

        plane = next(iter(observations[0].planes.values()))
        artifact = next(iter(plane.artifacts.values()))
        self.assertTrue(artifact.uri.startswith("cadc:TEST/_mock_calexp/"))

    def test_butler_uri_keyword_is_available(self):
        """{butler_uri} expands to empty string when there is no datastore."""
        from lsst.dax.obscore.caom_exporter import CaomExporter

        butler = self.make_populated_butler()
        config = self.make_caom_config()
        config.select_dataset_types(["_mock_calexp"])
        config.caom.dataset_types = {"_mock_calexp": config.caom.dataset_types["_mock_calexp"]}
        config.caom.artifact_uri_fmt = "cadc:TEST/{butler_uri}/{id}"

        observations = list(CaomExporter(butler, config).iter_observations())
        plane = next(iter(observations[0].planes.values()))
        artifact = next(iter(plane.artifacts.values()))
        self.assertTrue(artifact.uri.startswith("cadc:TEST//"))
```

- [ ] **Step 2: Run to verify they fail**

Run: `pytest tests/test_caom_exporter.py -q -k "uri"`

Expected: both FAIL — no warning is logged from `obscore_exporter`, and `{butler_uri}` raises `KeyError`.

- [ ] **Step 3: Add `_resolve_uris` to `ObscoreExporter`**

Add this method to `ObscoreExporter`, immediately after `iter_records`:

```python
    def _resolve_uris(self, refs: list[DatasetRef]) -> dict[Any, str]:
        """Resolve Butler URIs for a batch of dataset references.

        Parameters
        ----------
        refs : `list` [ `~lsst.daf.butler.DatasetRef` ]
            References to resolve.

        Returns
        -------
        uris : `dict` [ `~lsst.daf.butler.DatasetId`, `str` ]
            Mapping from dataset ID to primary URI. Datasets with no
            resolvable URI are absent, and the mapping is empty if the
            butler has no datastore.

        Notes
        -----
        Resolution is done in bulk rather than one call per dataset. A
        butler without a datastore, or files that are not present, produce
        a warning rather than an error, because the CAOM export can fall
        back to a configured URI template.
        """
        if not refs:
            return {}
        try:
            resolved = self.butler.get_many_uris(refs, predict=False)
        except Exception as exc:
            _LOG.warning("Could not resolve Butler URIs (%s); falling back to configured templates.", exc)
            return {}

        uris: dict[Any, str] = {}
        for ref, ref_uris in resolved.items():
            if ref_uris.primaryURI is not None:
                uris[ref.id] = str(ref_uris.primaryURI)
        return uris
```

- [ ] **Step 4: Use it in `CaomExporter`**

In `caom_exporter.py`, change `iter_observations` to resolve URIs in bulk and thread them through:

```python
    def iter_observations(self) -> Iterator[Any]:
        """Generate CAOM Observations for the configured dataset types.

        Yields
        ------
        observation : `caom2.Observation`
            One Observation per distinct expanded ``observation_id_fmt``.
        """
        observations: dict[str, Any] = {}
        state = _QueryState()
        pending: list[tuple[DatasetRef, Region | None, dict[str, Any]]] = list(
            self._obscore._iter_record_refs(state)
        )
        self._uris = self._obscore._resolve_uris([ref for ref, _, _ in pending])
        for ref, region, record in pending:
            self._add_record(observations, ref, region, record)
        yield from observations.values()
```

Add `self._uris: dict[Any, str] = {}` to `__init__` after `self._plane_owners`.

In `_format_keywords`, add the URI to the namespace before returning:

```python
        keywords.update(butler_uri=self._uris.get(ref.id, ""))
```

- [ ] **Step 5: Run the tests to verify they pass**

Run: `pytest tests/test_caom_exporter.py -q`

Expected: 7 passed.

- [ ] **Step 6: Confirm the ObsCore export is unaffected**

Run: `pytest tests/ -q`

Expected: all pass.

- [ ] **Step 7: Lint and commit**

Run: `ruff check python/lsst/dax/obscore/obscore_exporter.py python/lsst/dax/obscore/caom_exporter.py tests/test_caom_exporter.py && mypy python/lsst/dax/obscore/obscore_exporter.py python/lsst/dax/obscore/caom_exporter.py`

```bash
git add python/lsst/dax/obscore/obscore_exporter.py python/lsst/dax/obscore/caom_exporter.py tests/test_caom_exporter.py
git commit -m "Resolve Butler URIs in bulk for CAOM artifacts"
```

---

### Task 8: XML output, script layer and CLI command

**Files:**
- Modify: `python/lsst/dax/obscore/caom_exporter.py` (add `to_directory`)
- Create: `python/lsst/dax/obscore/script/obscore_export_caom.py`
- Modify: `python/lsst/dax/obscore/script/__init__.py`
- Modify: `python/lsst/dax/obscore/cli/cmd/commands.py`
- Test: `tests/test_caom_exporter.py`

**Interfaces:**
- Consumes: `CaomExporter.iter_observations` from Task 5.
- Produces: `CaomExporter.to_directory(destination: str, validate: bool = True) -> int` returning the number of documents written; `script.obscore_export_caom(repo, destination, config, where, collections, dataset_type)`; click command `export_caom`.

- [ ] **Step 1: Write the failing test**

Append to `CaomExporterTestCase`:

```python
    def test_to_directory_writes_valid_xml(self):
        """Documents are written per Observation and validate against XSD."""
        from lsst.dax.obscore.caom_exporter import CaomExporter

        butler = self.make_populated_butler()
        config = self.make_caom_config()
        destination = os.path.join(self.root, "caom")

        count = CaomExporter(butler, config).to_directory(destination)

        self.assertGreater(count, 0)
        written = sorted(os.listdir(destination))
        self.assertEqual(len(written), count)
        for name in written:
            self.assertTrue(name.endswith(".xml"))
            with open(os.path.join(destination, name), "rb") as handle:
                head = handle.read(512)
            self.assertIn(b"Observation", head)
```

- [ ] **Step 2: Run to verify it fails**

Run: `pytest tests/test_caom_exporter.py -q -k to_directory`

Expected: FAIL with `AttributeError: 'CaomExporter' object has no attribute 'to_directory'`.

- [ ] **Step 3: Add `to_directory`**

In `caom_exporter.py`, add `import os` and `import re` to the imports, and add this method immediately after `iter_observations`:

```python
    def to_directory(self, destination: str, validate: bool = True) -> int:
        """Write one CAOM XML document per Observation.

        Parameters
        ----------
        destination : `str`
            Directory to write into. Created if it does not exist.
        validate : `bool`, optional
            If `True`, validate each document against the bundled CAOM
            schema before writing.

        Returns
        -------
        count : `int`
            Number of documents written.
        """
        os.makedirs(destination, exist_ok=True)
        writer = caom2.ObservationWriter(validate=validate)

        count = 0
        for observation in self.iter_observations():
            filename = os.path.join(destination, f"{_safe_filename(observation.observation_id)}.xml")
            with open(filename, "wb") as handle:
                writer.write(observation, handle)
            count += 1

        _LOG.info("Wrote %d CAOM observation%s to %s", count, "" if count == 1 else "s", destination)
        return count
```

Add this module-level helper immediately before `class CaomExporter:`:

```python
def _safe_filename(observation_id: str) -> str:
    """Convert an observation identifier into a safe file name.

    Parameters
    ----------
    observation_id : `str`
        CAOM observation identifier.

    Returns
    -------
    name : `str`
        The identifier with characters that are awkward in file names
        replaced by underscores.
    """
    return re.sub(r"[^A-Za-z0-9._-]", "_", observation_id)
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `pytest tests/test_caom_exporter.py -q -k to_directory`

Expected: PASS. If schema validation rejects a document, read the error: it names the offending element, and the fix belongs in `_make_plane` or `_make_observation`, not in disabling validation.

The most likely failure is a `caom2.shape.Interval` for energy or time bounds needing an explicit `samples` list. If so, construct it as `caom2.shape.Interval(lower, upper, samples=[caom2.shape.SubInterval(lower, upper)])` in `_make_plane` and re-run.

- [ ] **Step 5: Add the script layer**

Create `python/lsst/dax/obscore/script/obscore_export_caom.py` with the GPL header, then:

```python
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
```

Add to `python/lsst/dax/obscore/script/__init__.py`, keeping alphabetical order (it goes after the `obscore_export` line):

```python
from .obscore_export_caom import obscore_export_caom
```

- [ ] **Step 6: Add the CLI command**

In `python/lsst/dax/obscore/cli/cmd/commands.py`, add this command immediately after the existing `export` command:

```python
@obscore.command(
    short_help="Export Butler datasets as CAOM observations",
    cls=ButlerCommand,
)
@repo_argument(required=True)
@destination_argument(
    required=True,
    help="DESTINATION is the directory to write the CAOM XML documents into.",
    type=MWPath(file_okay=False, dir_okay=True, writable=True),
)
@click.option(
    "--config",
    "-c",
    help="Location of the configuration file in YAML format, path or URL.",
    required=True,
)
@dataset_type_option(
    help=(
        "Comma-separated list of dataset types. "
        "If specified it must be a subset of the dataset types defined in the 'caom' "
        "section of the configuration file."
    )
)
@collections_option()
@where_option()
@options_file_option()
def export_caom(*args: Any, **kwargs: Any) -> None:
    """Export Butler datasets as CAOM observations, one XML document per
    observation.

    Requires the optional caom2 dependency, installed with
    'pip install lsst-dax-obscore[caom]'.
    """
    script.obscore_export_caom(*args, **kwargs)
```

Check `python/lsst/dax/obscore/cli/cmd/__init__.py` for an `__all__`; if the existing commands are listed there, add `"export_caom"` in sorted position.

- [ ] **Step 7: Verify the command is registered**

Run: `butler obscore --help`

Expected: `export-caom` appears in the subcommand list. Click converts the underscore in the function name to a hyphen.

Run: `butler obscore export-caom --help`

Expected: the options above are listed.

- [ ] **Step 8: Run the full suite, lint and commit**

Run: `pytest tests/ -q && ruff check python/ tests/ && ruff format --check python/ tests/ && mypy python/`

```bash
git add python/lsst/dax/obscore/caom_exporter.py python/lsst/dax/obscore/script/ \
    python/lsst/dax/obscore/cli/cmd/ tests/test_caom_exporter.py
git commit -m "Add butler obscore export-caom command"
```

---

### Task 9: Auxiliary artifacts

Attaches extra dataset types, such as `deep_coadd_n_image`, to the plane whose primary record shares their data ID. One extra query per auxiliary dataset type, not one lookup per plane.

**Files:**
- Modify: `python/lsst/dax/obscore/caom_exporter.py`
- Modify: `python/lsst/dax/obscore/tests.py` (`make_caom_config`)
- Test: `tests/test_caom_exporter.py`

**Interfaces:**
- Consumes: `CaomExporter._make_artifact` from Task 5.
- Produces: `CaomExporter._add_auxiliary_artifacts(observations, plane_data_ids) -> None`. `_add_record` records `plane_data_ids[(observation_id, product_id)] = ref.dataId`.

- [ ] **Step 1: Extend the test configuration**

In `make_caom_config` in `tests.py`, add an `auxiliary_datasets` entry to the `_mock_calexp` block:

```python
                        "auxiliary_datasets": {"_mock_src": "auxiliary"},
```

Check that `_mock_src` exists in `tests/data/hsc_gen3.yaml`:

Run: `grep -c "_mock_src" tests/data/hsc_gen3.yaml`

If the count is zero, pick a dataset type that is present and shares the visit and detector dimensions with `_mock_calexp`, and use that name instead throughout this task.

`test_calexp_observations` from Task 5 already asserts `assertGreaterEqual(len(plane.artifacts), 1)` rather than an exact count, so adding an auxiliary dataset type does not break it. Confirm that is still what the file says before continuing.

- [ ] **Step 2: Write the failing test**

Append to `CaomExporterTestCase`:

```python
    def test_auxiliary_artifacts(self):
        """Auxiliary dataset types add artifacts to the matching plane."""
        from lsst.dax.obscore.caom_exporter import CaomExporter

        butler = self.make_populated_butler()
        config = self.make_caom_config()
        config.select_dataset_types(["_mock_calexp"])
        config.caom.dataset_types = {"_mock_calexp": config.caom.dataset_types["_mock_calexp"]}

        observations = list(CaomExporter(butler, config).iter_observations())

        plane = next(iter(observations[0].planes.values()))
        product_types = {artifact.product_type for artifact in plane.artifacts.values()}
        self.assertIn(caom2.ProductType.THIS, product_types)
        self.assertIn(caom2.ProductType.AUXILIARY, product_types)
```

- [ ] **Step 3: Run to verify it fails**

Run: `pytest tests/test_caom_exporter.py -q -k auxiliary`

Expected: FAIL, because only the science artifact is present.

- [ ] **Step 4: Implement auxiliary artifact attachment**

In `_add_record`, after the artifact is added, record the plane's data ID. Change the signature to accept the mapping:

```python
    def _add_record(
        self,
        observations: dict[str, Any],
        plane_data_ids: dict[tuple[str, str], Any],
        ref: DatasetRef,
        region: Region | None,
        record: dict[str, Any],
    ) -> None:
```

and add, just before the artifact is created:

```python
        plane_data_ids[owner_key] = ref.dataId
```

Add this method to `CaomExporter`:

```python
    def _add_auxiliary_artifacts(
        self, observations: dict[str, Any], plane_data_ids: dict[tuple[str, str], Any]
    ) -> None:
        """Attach auxiliary dataset artifacts to the planes they belong to.

        Parameters
        ----------
        observations : `dict` [ `str`, `caom2.Observation` ]
            Observations built from the primary dataset types. Updated in
            place.
        plane_data_ids : `dict` [ `tuple` [ `str`, `str` ], \
                `~lsst.daf.butler.DataCoordinate` ]
            Data ID of the primary record for each plane, keyed by
            observation and product identifier.

        Notes
        -----
        Each auxiliary dataset type is queried once and its results indexed
        by data ID, rather than issuing a lookup for every plane.
        """
        collections = self.config.collections
        for dataset_type, config in self.caom_config.dataset_types.items():
            for auxiliary_type, product_type_name in config.auxiliary_datasets.items():
                by_data_id: dict[Any, DatasetRef] = {}
                try:
                    refs = self.butler.query_datasets(
                        auxiliary_type, collections=collections, explain=False
                    )
                except Exception as exc:
                    _LOG.warning("Could not query auxiliary dataset type %s: %s", auxiliary_type, exc)
                    continue
                for ref in refs:
                    by_data_id[ref.dataId] = ref

                uris = self._obscore._resolve_uris(list(by_data_id.values()))
                self._uris.update(uris)

                attached = 0
                for (observation_id, product_id), data_id in plane_data_ids.items():
                    if self._plane_owners.get((observation_id, product_id)) != dataset_type:
                        continue
                    ref = by_data_id.get(data_id)
                    if ref is None:
                        continue
                    keywords = self._format_keywords(ref, {})
                    artifact = self._make_artifact(config, keywords, product_type_name)
                    if artifact is not None:
                        observations[observation_id].planes[product_id].artifacts[artifact.uri] = artifact
                        attached += 1

                _LOG.info(
                    "Attached %d auxiliary artifact%s of type %s",
                    attached,
                    "" if attached == 1 else "s",
                    auxiliary_type,
                )
```

Update `iter_observations` to thread the mapping through and call the new method:

```python
        observations: dict[str, Any] = {}
        plane_data_ids: dict[tuple[str, str], Any] = {}
        state = _QueryState()
        pending: list[tuple[DatasetRef, Region | None, dict[str, Any]]] = list(
            self._obscore._iter_record_refs(state)
        )
        self._uris = self._obscore._resolve_uris([ref for ref, _, _ in pending])
        for ref, region, record in pending:
            self._add_record(observations, plane_data_ids, ref, region, record)
        self._add_auxiliary_artifacts(observations, plane_data_ids)
        yield from observations.values()
```

Note that `_format_keywords` is called with an empty record for auxiliary datasets, because there is no ObsCore record for them; `artifact_uri_fmt` templates that reference ObsCore columns will fail for auxiliary types. Document that in the configuration docstring in Task 10.

- [ ] **Step 5: Run the tests to verify they pass**

Run: `pytest tests/test_caom_exporter.py -q`

Expected: 9 passed.

- [ ] **Step 6: Lint and commit**

Run: `ruff check python/lsst/dax/obscore/caom_exporter.py python/lsst/dax/obscore/tests.py tests/test_caom_exporter.py && mypy python/lsst/dax/obscore/caom_exporter.py`

```bash
git add python/lsst/dax/obscore/caom_exporter.py python/lsst/dax/obscore/tests.py tests/test_caom_exporter.py
git commit -m "Attach auxiliary dataset artifacts to CAOM planes"
```

---

### Task 10: DP1 configuration and documentation

**Files:**
- Modify: `configs/dp1.yaml`
- Create: `doc/lsst.dax.obscore/caom-export.rst`
- Modify: `doc/lsst.dax.obscore/index.rst`
- Create: `doc/changes/DM-55951.feature.rst`

**Interfaces:**
- Consumes: everything above.
- Produces: no code.

- [ ] **Step 1: Add the `caom` block to `configs/dp1.yaml`**

Append at the end of the file, after the existing `origin:` block:

```yaml
caom:
  telescope_name: "8.4-meter Simonyi Survey Telescope"
  proposal_id: LSST
  em_band: OPTICAL
  provenance:
    name: "LSST Science Pipelines"
    version: "v29.0"
    project: DRP
  read_groups:
    - "ivo://cadc.nrc.ca/gms?LSST_CADC"
  intent: science
  observation_type: science
  content_type: application/fits
  artifact_uri_fmt: "cadc:LSST/{dataset_type}/{id}"
  dataset_types:
    raw:
      observation_id_fmt: "{records[exposure].obs_id}"
      product_id_fmt: "raw-{records[detector].full_name}"
      observation_type_fmt: "{records[exposure].observation_type}"
      s_pixel_scale: 0.2
    visit_image:
      observation_id_fmt: "{records[visit].name}"
      product_id_fmt: "visit_image-{records[detector].full_name}"
      s_pixel_scale: 0.2
    difference_image:
      observation_id_fmt: "{records[visit].name}"
      product_id_fmt: "difference_image-{records[detector].full_name}"
      s_pixel_scale: 0.2
    deep_coadd:
      observation_id_fmt: "{skymap}-{tract}"
      product_id_fmt: "deep_coadd-{patch}-{band}"
      derived: true
      algorithm: "lsst.deep_coadd.DRP.DP1.DM-51335"
      s_pixel_scale: 0.2
    template_coadd:
      observation_id_fmt: "{skymap}-{tract}"
      product_id_fmt: "template_coadd-{patch}-{band}"
      derived: true
      algorithm: "lsst.template_coadd.DRP.DP1.DM-51335"
      s_pixel_scale: 0.2
```

- [ ] **Step 2: Verify the configuration validates**

Run:

```bash
python -c "
from lsst.daf.butler import Config
from lsst.dax.obscore import ExporterConfig
cfg = ExporterConfig.model_validate(Config('configs/dp1.yaml'))
assert cfg.caom is not None
print(sorted(cfg.caom.dataset_types))
"
```

Expected: the five dataset type names, in sorted order.

Run: `yamllint configs/dp1.yaml`

Expected: clean.

- [ ] **Step 3: Write the documentation page**

Create `doc/lsst.dax.obscore/caom-export.rst`. Use one sentence per line, American English spelling. It must cover, each as its own section:

1. What the command does and how to run it, showing a complete `butler obscore export-caom` invocation against `configs/dp1.yaml`.
2. The optional `caom2` dependency and the `lsst-dax-obscore[caom]` extra.
3. Every key in the `caom` block, with its meaning and default, taken from the docstrings in `caom_config.py`.
4. The table of CAOM elements derived from existing ObsCore configuration, copied from the spec's "Configuration" section.
5. **Observation merging semantics**, stated explicitly: dataset types whose `observation_id_fmt` templates expand to the same value contribute planes to a single Observation; the first dataset type to reach an observation ID sets the observation-level attributes and a later disagreement is warned about and ignored; a duplicate product ID between dataset types is an error.
6. The templating namespace: the data ID mapping, `records`, `id`, `run`, `dataset_type`, `butler_uri`, and every ObsCore column. Note that auxiliary dataset types have no ObsCore record, so their `artifact_uri_fmt` may only use the data ID, `id`, `run`, `dataset_type` and `butler_uri`.
7. Known gaps, copied from the spec: the dataset DOI has no CAOM home, `contentLength` is unset, there are no `Part` or `Chunk` elements, and `s_pixel_scale` is a placeholder that will move to the ObsCore `DatasetTypeConfig`.
8. The memory characteristic: all observations accumulate before writing, so use `--where` to chunk a large export.

Add `caom-export` to the toctree in `doc/lsst.dax.obscore/index.rst`.

- [ ] **Step 4: Add a change log fragment**

Check the naming convention:

Run: `ls doc/changes/`

Create `doc/changes/DM-55951.feature.rst` following whatever convention the existing files use, with content:

```
Added a new ``butler obscore export-caom`` command that exports Butler datasets as CAOM observations.
The command is driven by a new ``caom`` section in the existing ObsCore configuration files, so the CAOM output is built from the same records as the ObsCore output.
It requires the optional ``caom2`` dependency, installed with ``pip install lsst-dax-obscore[caom]``.
```

- [ ] **Step 5: Build the docs**

Run: `package-docs build` from the package root, or `pytest tests/ -q` if the Sphinx toolchain is not set up locally.

Expected: no Sphinx errors referencing the new page.

- [ ] **Step 6: Full verification**

Run: `pytest tests/ -q && ruff check python/ tests/ && ruff format --check python/ tests/ && mypy python/ && yamllint configs/`

Expected: all clean.

- [ ] **Step 7: Commit**

```bash
git add configs/dp1.yaml doc/
git commit -m "Add DP1 caom configuration and CAOM export documentation"
```

---

## Follow-on work

Not planned in detail here, because phase 4 may change the answers.

**Phase 4, CADC review.** Produce output from a real DP1 repository and send it to CADC. Settle: whether tract-level observations with patch and band planes is the grouping they want; the correct `Artifact.uri` namespace; whether the dataset DOI needs a home; and whether `contentLength` is required at ingest.

**Phase 5, move `s_pixel_scale` to ObsCore.** File a `daf_butler` ticket adding `s_pixel_scale` to `DatasetTypeConfig` beside `s_xel`, and surfacing it on the ObsCore record. Then change `CaomConfig.pixel_scale_for` to read the record instead of the placeholder key, remove `s_pixel_scale` from `CaomDatasetTypeConfig`, and move the values in `configs/dp1.yaml` from the `caom` block to the top-level `dataset_types` block. The single accessor is why this is a small change.

**Unrelated cleanup, worth its own ticket.** `use_butler_uri` appears in nine files under `configs/` and in `python/lsst/dax/obscore/tests.py`, and no code reads it. Pydantic silently discards it. Either implement it or remove it.
