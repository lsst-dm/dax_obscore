# `butler obscore export-regions` design

Ticket: DM-53569

## Summary

Add a new `obscore export-regions` subcommand that walks a butler collection,
selects datasets via a `--where` clause, extracts the spatial region for each
dataset, deduplicates, converts each region to an IVOA STC-S string via
`lsst.sphgeom.Region.to_ivoa_stcs()`, and writes a FITS binary table whose
single `s_region` column is in the form expected as input to the TOPCAT STILTS
`mocshape` command.

## Scope

In scope:

- A new click subcommand `export-regions` under the existing `obscore` group.
- A script function `obscore_export_regions(...)` that does the work.
- Support for dataset types whose `DimensionGroup.region_dimension` resolves
  directly — i.e., `visit`-, `visit_detector_region`-, `tract`-, and
  `patch`-style spatial dimensions.
- Tests covering the script function and a CLI smoke test.

Out of scope:

- Computing the MOC itself. The output is the *input* to `stilts mocshape`,
  not a MOC FITS file.
- ObsCore configuration. This command does not consume an `ExporterConfig`
  YAML and does not reuse `SIAv2Handler` or `ObscoreExporter`.
- Exposure-only dataset types whose region must be synthesized via a fallback
  to `visit` or `visit_detector_region`. These will not be exercised by users
  of this command and are not supported.
- `Box` regions. The butler does not produce them for the supported dataset
  types; if `to_ivoa_stcs()` ever raises, the command fails hard.

## CLI surface

```
butler obscore export-regions REPO DESTINATION \
    --dataset-type <name> \
    --collections <coll>[,<coll>...] \
    --where "<expression>"
```

Implemented in `python/lsst/dax/obscore/cli/cmd/commands.py` as a new
`@obscore.command` decorated with `cls=ButlerCommand`. Uses these standard
daf_butler click options: `repo_argument`, `destination_argument`,
`collections_option`, `where_option`, `options_file_option`. The
`--dataset-type` option is declared inline as a single required string
(`type=str`, no `multiple=True`, no `split_commas`) — this differs from the
existing `export`/`siav2` subcommands, which accept multiple comma-separated
dataset types.

`DESTINATION` is the path to the output FITS file.

## Architecture

Two new pieces, mirroring the file layout used by the existing subcommands:

1. **`python/lsst/dax/obscore/script/obscore_export_regions.py`** — module
   containing the `obscore_export_regions(repo, destination, dataset_type,
   collections, where)` function. This is where all real work happens.

2. **A new subcommand in `python/lsst/dax/obscore/cli/cmd/commands.py`** —
   thin click wrapper that calls the script function.

No new public APIs in `lsst.dax.obscore` outside the script module.

## Region selection

The dataset type's `DimensionGroup.region_dimension` is the single source of
truth for which dimension carries the region:

```python
dataset_type = butler.get_dataset_type(name)
region_dim = dataset_type.dimensions.region_dimension
if region_dim is None:
    raise ValueError(
        f"Dataset type {name!r} has no spatial region dimension."
    )
region_key = f"{region_dim}.region"
```

This naturally yields `visit_detector_region` for `(visit, detector)` dataset
types, `visit` for `(visit,)` dataset types, `patch` for `(skymap, tract,
patch)` dataset types, and `tract` for `(skymap, tract)` dataset types.

We do **not** replicate the SIAv2 fallback that maps an `exposure`-only
dataset type to `visit_detector_region`, because users of this command will
only run it on visit- or patch-bearing dataset types.

## Querying and deduplication

The query follows the established `obscore_exporter` pattern:

```python
with butler.query() as query:
    q = query.join_dataset_search(dataset_type.name, collections=collections)
    if where:
        q = q.where(where)
    result = q.general(
        dataset_type.dimensions,
        region_key,
        dataset_fields={dataset_type.name: ...},
        find_first=True,
    )
    required_keys = sorted(
        butler.dimensions[region_dim].minimal_group.required
    )
    seen: set[tuple] = set()
    regions: list[Region] = []
    for dataId, _refs, raw_row in result.iter_tuples(dataset_type):
        key = tuple(dataId[k] for k in required_keys)
        if key in seen:
            continue
        seen.add(key)
        regions.append(raw_row[region_key])
```

Notes:

- Dedup key is the dataId restricted to the required keys of the region
  dimension (e.g., `(instrument, visit, detector)` for `visit_detector_region`,
  `(skymap, tract, patch)` for `patch`). The required-key list is sorted
  before iteration so the tuple ordering is stable across runs.
- The exact API used to extract those required keys is the
  `DimensionElement.minimal_group.required` set; minor adjustments are
  acceptable during implementation as long as the dedup semantics match.
- We accumulate regions in a Python list. With 128 GB hosts and at most a
  few million regions, streaming is unnecessary.

## STC-S conversion

For each unique region:

```python
stcs_strings.append(region.to_ivoa_stcs())  # default frame=ICRS
```

If `to_ivoa_stcs()` raises (e.g., for a `Box`), the exception propagates and
the command fails. This is the agreed behavior.

## FITS output

A single-extension binary table FITS file written with `astropy.io.fits`.

- Extension name: `EXTNAME = 'REGIONS'`.
- One column:
  - Name: `s_region`
  - Type: string. Implementation may use either an astropy variable-length
    character column (`PA(maxlen)`) or a fixed-width column (`An`) sized to
    the longest STC-S string produced; choose whichever astropy and STILTS
    handle most cleanly during implementation.
- Column-level header keywords on the BinTableHDU header (column index 1):
  - `TUTYP1 = 'obscore:Char.SpatialAxis.Coverage.Support.Area'`
  - `TUCD1  = 'pos.outline;obs.field'`
  - `TXTYP1 = 'stc-s'`
- A `HISTORY` line recording the dataset type, collections, and where clause
  used to produce the file (analogous to `query_url` for `siav2`).

The `TXTYP1 = 'stc-s'` keyword is what allows STILTS `mocshape` to detect the
column as STC-S without an explicit `shape=stc-s` argument; the `TUTYP1` and
`TUCD1` values match ObsCore convention for `s_region`.

Empty input (no datasets match) produces a valid FITS file with zero rows.

Output is written via `hdulist.writeto(destination, overwrite=True)`.

## Downstream usage

```
stilts mocshape in=regions.fits coords=s_region order=10 out=coverage.moc.fits
```

STILTS auto-detects the STC-S column from `TXTYP1`.

## Testing

New test module `tests/test_export_regions.py`, reusing the existing daf_butler
test repository fixture used by the other obscore tests.

Cases:

1. Visit/detector dataset type — happy path. Writes a FITS file; row count
   equals the number of distinct `(instrument, visit, detector)` triples;
   every `s_region` cell parses as STC-S; header has `TXTYP1='stc-s'`,
   `TUTYP1='obscore:Char.SpatialAxis.Coverage.Support.Area'`,
   `TUCD1='pos.outline;obs.field'`.
2. Patch-bearing dataset type — happy path. Returns distinct patch regions.
3. Deduplication. N datasets sharing K distinct region-dim keys produce K
   rows.
4. `--where` filtering. A subsetting where clause produces the expected
   reduced row count.
5. Error: dataset type with no spatial dimension raises `ValueError`.
6. Empty result: no datasets match the where clause → valid FITS file with
   zero rows.

CLI smoke test: invoke `obscore export-regions` through Click's `CliRunner`
against the same fixture and assert exit code 0 and output file presence.

We do **not** re-test sphgeom's STC-S correctness; that lives on the
`tickets/DM-53569` branch of `sphgeom`.

## Dependencies

- `lsst.sphgeom` on the `tickets/DM-53569` branch (provides
  `Region.to_ivoa_stcs()`).
- `lsst.daf.butler` — already a runtime dependency.
- `astropy.io.fits` — already a transitive dependency via existing modules.

No new third-party dependencies.
