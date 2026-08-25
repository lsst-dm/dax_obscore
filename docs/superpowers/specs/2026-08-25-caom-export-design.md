# CAOM export from Butler (`butler obscore export-caom`)

Ticket: DM-55951

Date: 2026-08-25

## Background

CADC maintain a standalone Butler-to-CAOM export script at
[`lsst2caom2/butler2caom.py`](https://github.com/ijiraq/lsst2caom2/blob/main/lsst2caom2/butler2caom.py).
That script builds a CAOM `Observation` for a single LSSTComCam DP1 `deep_coadd` and prints it.
It is 333 lines, of which a large fraction is module-level constants describing DP1:
telescope geolocation, instrument name, proposal, bandpass throughput ranges, pixel scale, coadd dimensions, and the CAOM collection name.
It reads each image with `butler.get()` in order to recover the WCS, and issues an authenticated HTTP `HEAD` per file to obtain the content length.

dax_obscore already computes almost every quantity that script derives, and it does so from registry records alone.
`s_region`, `s_ra`, `s_dec`, `s_fov`, `t_min`, `t_max`, `t_exptime`, `em_min`, `em_max`, `target_name`, `instrument_name`, `calib_level` and `dataproduct_type` are all produced by `lsst.daf.butler.registry.obscore.RecordFactory` without touching a file.
Porting the CADC script into dax_obscore therefore removes the file reads, removes the hardcoded DP1 constants, and — most importantly — makes the CAOM output a second rendering of the same records that feed the ObsCore output, rather than an independent description that can drift away from it.

## Goals

Add a `butler obscore export-caom` subcommand that writes CAOM XML documents from a Butler repository, driven by the existing ObsCore configuration files extended with a CAOM section.

Specifically:

- No file access. Everything derives from registry records and configuration.
- No divergence. The CAOM output is built from the same ObsCore records the `obscore export` command emits, produced by the same query.
- Generic over dataset types. Any dataset type given a CAOM entry in the configuration can be exported; there is no DP1-specific or coadd-specific code.
- Fix gaps in the ObsCore layer rather than working around them in the CAOM layer. Any quantity that is an ObsCore concept is populated on the ObsCore side and consumed from the record, so that the ObsCore export benefits too. CAOM-only configuration is reserved for genuinely CAOM-specific concepts: the proposal, the read groups, the algorithm, and the observation and plane grouping templates.

## Non-goals

This work is a prototype intended to trigger discussion with CADC, and the decisions recorded here are expected to be revised.
The following are deliberately out of scope.

- `Part` and `Chunk` elements. CAOM 2.5 removes them, and the discovery metadata CADC needs sits at `Plane` level.
- Uploading to a `caom2repo` service. The command writes files.
- `Artifact.contentLength`. See "Known gaps" below.

## Design decisions

| Decision | Choice |
| --- | --- |
| Scope | Generic over dataset types configured with a CAOM entry |
| Output | A directory of CAOM XML documents, one per `Observation` |
| Observation identity | Config-driven template, shared across dataset types |
| Cross-dataset-type merging | Yes — planes from several dataset types merge into one `Observation` |
| CAOM depth | `Observation` → `Plane` → `Artifact`; no `Part`/`Chunk` |
| Config location | A single `caom` property on `ExporterConfig`, in the existing ObsCore config file |
| `caom2` dependency | Optional extra, imported lazily |
| Auxiliary artifacts | In scope, listed per dataset type in configuration |
| `Artifact.contentLength` | Unset in the prototype |
| ObsCore gaps | Fixed on the ObsCore side and read from the record, never worked around in the CAOM layer |

## Architecture

Three new modules, following the flat-module layout the package already uses.

| File | Contents | Imports `caom2` |
| --- | --- | --- |
| `python/lsst/dax/obscore/caom_config.py` | Pydantic models for the `caom` configuration block | No |
| `python/lsst/dax/obscore/caom_exporter.py` | `CaomExporter`, which builds and writes the documents | Yes, lazily |
| `python/lsst/dax/obscore/script/obscore_export_caom.py` | Script layer invoked by the CLI | No |

`ExporterConfig` gains a `caom: CaomConfig | None` field.
`ExporterConfig` is imported by the RSP SIA service, so `caom_config.py` must not import `caom2` under any circumstance.
Only `caom_exporter.py` imports it, inside a `try`/`except ImportError` that records the failure and raises a message naming the `lsst-dax-obscore[caom]` extra when the exporter is actually used.

### Reuse of `ObscoreExporter`

`CaomExporter` does not reimplement the registry query.
It consumes records produced by `ObscoreExporter`, which requires one contained refactor of `obscore_exporter.py`.

The query and record-construction loop currently lives inside `ObscoreExporter._make_record_batches`, which interleaves it with pyarrow batching and with overflow accounting for the SIAv2 `QUERY_STATUS` report.
That loop is extracted into a generator yielding one tuple per matching dataset:

```python
def _iter_record_refs(
    self, limit: int | None = None
) -> Iterator[tuple[DatasetRef, Region | None, ResourcePath | None, dict[str, Any]]]:
```

`_make_record_batches` becomes a thin wrapper that batches the generator's output.
Overflow state is carried on a small mutable holder passed into the generator so that `to_votable`'s existing `OVERFLOW` reporting is unchanged.
`iter_records` and the three `to_*` methods keep their current signatures and behavior, and the existing tests in `tests/test_exporter.py` cover the refactor.

Yielding `region` rather than re-parsing the `s_region` string matters: the loop already holds the `lsst.sphgeom.Region`, and a sphgeom `ConvexPolygon` converts directly into `caom2.shape.Polygon` vertices with no STC-S round-trip.
Yielding `ref` gives `CaomExporter` the same format vocabulary that `obs_id_fmt` already uses, and the resolved `uri` supplies `Artifact.uri` as described under "ObsCore-side fixes".

## Configuration

A single `caom` property on `ExporterConfig`, held in the existing ObsCore configuration file.
Only dataset types listed under `caom.dataset_types` are exported; the exporter narrows the ObsCore configuration to that set before querying.
Keys under `caom.dataset_types` must be a subset of the top-level `dataset_types` keys, and validation fails otherwise.

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
  dataset_types:
    raw:
      observation_id_fmt: "{records[exposure].obs_id}"
      product_id_fmt: "raw-{records[detector].full_name}"
      observation_type_fmt: "{records[exposure].observation_type}"
    visit_image:
      observation_id_fmt: "{records[visit].name}"
      product_id_fmt: "visit_image-{records[detector].full_name}"
    difference_image:
      observation_id_fmt: "{records[visit].name}"
      product_id_fmt: "difference_image-{records[detector].full_name}"
    deep_coadd:
      observation_id_fmt: "{skymap}-{tract}"
      product_id_fmt: "deep_coadd-{patch}-{band}"
      algorithm: "lsst.deep_coadd.DRP.DP1.DM-51335"
      derived: true
      auxiliary_datasets:
        deep_coadd_n_image: auxiliary
        deep_coadd_background: auxiliary
```

Every value that the ObsCore configuration already carries is derived rather than restated.

| CAOM element | Source |
| --- | --- |
| `Observation.collection` | `obs_collection` |
| `Observation.instrument` | record `instrument_name`, which honors `fallback_instrument` |
| `Observation.target` | record `target_name` |
| `Observation.metaRelease` | `origin.publication_date` |
| `Telescope.name` | `caom.telescope_name`, defaulting to `facility_name` |
| `Telescope.geoLocation{X,Y,Z}` | `EarthLocation.of_site()` on the facility name |
| `Proposal.title` | `origin.title` |
| `Proposal.project` | `obs_collection` |
| `Provenance.producer` | `origin.publisher` |
| `Provenance.reference` | `origin.reference_url` |
| `Provenance.runID` | `ref.run` |
| `Plane.calibrationLevel` | `dataset_types.<X>.calib_level` |
| `Plane.dataProductType` | `dataset_types.<X>.dataproduct_type` |
| `Plane.observable` | `dataset_types.<X>.o_ucd` |
| `Position.bounds` | the `lsst.sphgeom.Region` from the query |
| `Position.dimension` | `dataset_types.<X>.s_xel` |
| `Position.sampleSize` | record `s_pixel_scale`, converted from arcsec to degrees |
| `Position.resolution` | record `s_resolution`, when populated |
| `Energy.bounds` | `spectral_ranges`, via record `em_min`/`em_max` |
| `Energy.bandpassName` | record `em_filter_name` |
| `Time.bounds`, `Time.exposure` | record `t_min`/`t_max`/`t_exptime` |
| `Artifact.uri` | the Butler URI, optionally rewritten by `caom.artifact_uri_fmt` |

This is where `BANDPASS_TOTAL_THROUGHPUT`, `PIXEL_SCALE`, `COADD_DIMENSION`, `TELESCOPE`, `INSTRUMENT`, `OBSTYPE`, `CAOM_COLLECTION` and the remaining upstream module constants dissolve into configuration that already exists.

### Notes on individual fields

`EarthLocation.of_site` resolves `"Rubin:Simonyi"`, so the facility name is sufficient and no separate site key is needed.
Astropy does not know every AAS facility code, so an explicit `caom.geo_location` override is available for facilities it cannot resolve.

`Provenance.version` is the version of the software that produced the data, a fixed release-level fact such as `v29.0` for DP1 and `v30.0` for DP2.
It is not a template.
The per-execution value belongs in `Provenance.runID`, which is set to `ref.run`.
The upstream script derived both from positional slices of the run string, taking `UPath(ref.run).parts[-1]` as the version and `parts[2]` as the project; those slices happen to produce sensible values for DP1 run names and produce nonsense for anything else, so they are replaced by explicit configuration.

`Observation.type` is the `OBSTYPE`-style description of the exposure.
The `exposure` dimension carries `observation_type` with values such as `science`, `bias` and `dark`, so exposure-based dataset types read it from the record through `observation_type_fmt`.
The `visit` dimension has no such field because a visit is always an on-sky science observation, so visit-based and coadd dataset types take the literal default, which is `science` rather than the upstream `OBJECT` so that the vocabulary matches the Rubin data model across all dataset types.
`intent` follows the same pattern, since a `raw` bias is a calibration rather than a science observation.

`content_type` varies by dataset type and is not always FITS; the DP2 DOI configuration records `application/vnd.apache.parquet` for catalog products alongside `application/fits` for images.
It is therefore a per-dataset-type key with a top-level default.
It cannot reuse ObsCore's `access_format`, which in `configs/dp1.yaml` describes the DataLink response rather than the file.

`Position.sampleSize` is the pixel scale and `Position.resolution` is the spatial resolution; they are different quantities and are handled separately in "ObsCore-side fixes" below.
Neither may be derived from `s_fov`, which is a bounding-circle diameter rather than an image width.

`algorithm` applies only when `derived` is true.
The `caom2` library forces `SimpleObservation` to use the algorithm name `exposure` and requires `DerivedObservation` to use anything else, so a configured `algorithm` on a non-derived dataset type is a validation error.

## ObsCore-side fixes

Three quantities the CAOM export needs are ObsCore concepts that the ObsCore layer does not currently supply.
Each is addressed on the ObsCore side, so that whatever is fixed benefits the ObsCore export as well and nothing is reintroduced as a CAOM-only workaround.
Two are fixed here; the third is recorded as a gap that registry records cannot close.

**Pixel scale.**
CAOM `Position.sampleSize` is the pixel scale.
ObsCore has no standard column for it, and it must not be conflated with `s_resolution`, which is the PSF FWHM.
A new per-dataset-type key `s_pixel_scale`, in arcsec, is added to `DatasetTypeConfig` in `lsst.daf.butler.registry.obscore`, mirroring the existing `s_xel`, and surfaced on the record.
This needs a companion `daf_butler` ticket.
It sits with `s_xel` in the ObsCore configuration rather than in the `caom` block, so a single statement of the dataset type's geometry serves both exports.

```yaml
dataset_types:
  deep_coadd:
    s_xel: [3400, 3400]
    s_pixel_scale: 0.2      # arcsec, new
```

**`s_resolution`.**
This is a standard ObsCore column that `RecordFactory` never populates, so DP1's published ObsCore output has it `NULL` today.
It is the PSF FWHM, which is measured per dataset and is not available from registry records, so it cannot be filled from configuration.
It stays `NULL`, and CAOM leaves `Position.resolution` unset.
This is recorded here as a known ObsCore gap rather than papered over with a pixel scale on the CAOM side, which would be the wrong quantity.

**The Butler file URI.**
`Artifact.uri` needs the location of the file itself.
ObsCore `access_url` is a DataLink URL whenever `use_butler_uri` is false, as it is in `configs/dp1.yaml`, so it cannot serve.
The URI is resolved in bulk from the Butler datastore records for the refs returned by each query, not by a call per dataset, and is exposed as a fourth element of the shared generator's tuple so both exports can use it.
This is the same access path the eventual `file_size` work described under "Known gaps" needs.
`caom.artifact_uri_fmt` then rewrites that URI into the namespace CADC ingest expects, which is a genuinely CAOM-specific concern and correctly belongs in the `caom` block.

`content_type` is a fourth candidate.
It is not a standard ObsCore column, and `access_format` describes the DataLink response rather than the file, so it stays in the `caom` block for now.
Whether it should instead join `s_xel` and `s_pixel_scale` as a per-dataset-type ObsCore key is an open question for review.

## Data flow

1. Validate that `caom.dataset_types` is a subset of `dataset_types`, then narrow the ObsCore configuration to those dataset types.
2. Iterate the shared generator once per dataset type. For each `(ref, region, uri, record)`:
   - Build a single format namespace containing the dataId mapping, `records`, `id`, `run`, `dataset_type`, and every finished ObsCore column, so that `{tract}`, `{records[visit].name}`, `{obs_id}` and `{lsst_patch}` are all usable in any template.
   - Expand `observation_id_fmt` and `product_id_fmt`.
   - Get or create the `Observation` in an in-memory mapping keyed by observation ID, then get or create the `Plane` within it.
   - Build the `Artifact` from the ref and append it to the plane.
   - Fill `Plane.position`, `Plane.energy` and `Plane.time` from `region` and the record columns when the plane is first created.
3. For each auxiliary dataset type, run one additional query using the same collections and `where` constraint, index the results by dataId, and attach the resulting artifacts to the plane whose primary record shares that dataId. This is one query per auxiliary dataset type, not a lookup per plane.
4. Write one `<observation_id>.xml` per `Observation` into the destination directory using `caom2.ObservationWriter`.

### Observation merging semantics

Records from different dataset types that expand to the same observation ID merge into a single `Observation` whose planes sit at different calibration levels.
For a visit, `raw` contributes a calibration-level-1 plane, `visit_image` a level-2 plane and `difference_image` a level-3 plane, and CADC receives one ingestable document per visit.
This is the canonical CAOM shape and is the reason the exporter accumulates rather than streams.

The grouping is entirely a property of the templates.
For DP1 the intended shapes are one `Observation` per visit with one plane per detector for visit-based data, and one `Observation` per tract with one plane per patch and band for coadds.
Because a tract observation ID and a visit observation ID can never collide, coadds never merge with visit data.

Two conflict rules apply, and both are logged.

- Dataset types sharing an observation ID but disagreeing on instrument, intent, type or `derived` produce a warning, and the first record encountered wins.
- Two records producing the same `(observation_id, product_id)` pair from different dataset types are a configuration error and abort the export, since distinct dataset types are expected to carry distinct product IDs.

### Region conversion

`Plane.position.bounds` is converted from the `lsst.sphgeom.Region` directly.
`ConvexPolygon` becomes a `caom2.shape.Polygon`, `Circle` becomes a `caom2.shape.Circle` and `Box` becomes a `caom2.shape.Box`.
`caom2.shape.Polygon` requires both a vertex list and a `MultiPolygon` sample set, and validates vertex ordering, so the converter constructs both.
Any other region type, including union regions, produces a warning and a plane with no position rather than an aborted export.

### Memory

All observations accumulate before anything is written, because a single observation can gather planes from several dataset types that are queried in separate passes.
At DP1 scale this is a few thousand observations and is not a concern.
It is the known scaling limit of the prototype, and `--where` is the escape hatch for chunking a larger export.
A future version can sort by observation ID and stream.

## Command-line interface

```
butler obscore export-caom REPO DESTINATION -c CONFIG
    [--dataset-type ...] [--collections ...] [--where ...]
```

`DESTINATION` is a directory rather than a file, declared as `MWPath(file_okay=False, dir_okay=True, writable=True)` and created when absent.
The remaining options mirror `obscore export` exactly, so both commands accept the same configuration file and the same filters.
`--dataset-type` narrows further within `caom.dataset_types` and reports an error when given a dataset type that has no CAOM entry.

## Testing

`tests/test_caom_config.py` covers configuration validation and does not require `caom2`.
It checks that a `caom.dataset_types` key absent from `dataset_types` is rejected, that `algorithm` on a non-derived dataset type is rejected, and that the defaults resolved from `origin`, `facility_name` and `obs_collection` are correct.

`tests/test_caom_exporter.py` extends `DaxObsCoreTestMixin` with a `caom` block covering the existing `_mock_calexp` and `_mock_deepCoadd` mock dataset types.
It asserts observation grouping, plane merging across dataset types, the position, energy and time values on a plane, and the artifact list.
The module is guarded with `unittest.skipUnless` on `caom2` being importable.

`tests/test_exporter.py` is unchanged and covers the `_make_record_batches` refactor, which is expected to produce identical output.

## Packaging

`caom2` is an optional dependency.

- `pyproject.toml` gains `[project.optional-dependencies] caom = ["caom2"]`.
- `requirements-caom.txt` pins it for CI, and the build workflow installs it so the exporter tests run.
- It is not added to `ups/dax_obscore.table`, so EUPS stack users without it receive the lazy-import error message rather than a broken build.

## Documentation

A page under `doc/lsst.dax.obscore/` documents the `caom` configuration block, the field-derivation table above, and the observation-merging semantics.
`configs/dp1.yaml` gains a `caom` block as the worked example.

## Known gaps and open items for CADC

These are recorded rather than solved, and are the intended subject of the discussion this prototype exists to start.

- **The dataset DOI has no home in CAOM.** CAOM 2.4 has no DOI field. `Plane.creatorID` is an IVOA dataset identifier and is not a DOI. `origin.citation` is therefore not exported.
- **`Artifact.contentLength` is unset.** The Butler file datastore records carry `file_size`, which is the fastest source and requires no network access, but that value is sometimes `-1`. The intended resolution is to read the datastore records and fall back to `ResourcePath.size()` for the `-1` cases behind an opt-in flag, since the fallback issues one HTTP `HEAD` per file and needs credentials for the target store.
- **No `Part` or `Chunk` elements.** CAOM 2.5 removes them. A per-dataset WCS could be reconstructed from the Butler region if they turn out to be needed before then, still without reading a file.
- **The grouping is a proposal.** Whether tract-level observations with patch and band planes, and visit-level observations with detector planes, are the shapes CADC wants is the main question to settle.
- **ObsCore `s_resolution` is unpopulated.** See "ObsCore-side fixes". Filling it requires PSF FWHM per dataset, which the registry does not hold.
- **`Artifact.uri` rewriting.** The upstream script stripped a `/raven/files/` prefix from the Butler URI path. `caom.artifact_uri_fmt` generalizes that, but the correct target URI form for CADC ingest needs confirming.
