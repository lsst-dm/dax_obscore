.. _lsst.dax.obscore-caom-export:

###########
CAOM export
###########

``butler obscore export-caom`` exports Butler datasets as `Common Archive Observation Model <https://www.opencadc.org/caom2/>`_ (CAOM) observations, writing one XML document per observation.

The command reads nothing but registry records and configuration.
It never opens a data file.
Every value it writes is either taken from the ObsCore record that ``butler obscore export`` would produce for the same dataset, or from the ``caom`` section of the same configuration file.
Building both exports from one set of records is deliberate: it prevents the CAOM description of a dataset from drifting away from the ObsCore description of it.

Running the command
===================

.. code-block:: bash

   butler obscore export-caom /repo/dp1 caom-output -c configs/dp1.yaml

``DESTINATION`` is a directory rather than a file, and is created if it does not exist.
Each observation is written to ``<observation_id>.xml`` within it, with characters that are awkward in a file name replaced by underscores.

The ``--dataset-type``, ``--collections`` and ``--where`` options behave as they do for ``butler obscore export``, so both commands accept the same configuration file and the same filters.
``--dataset-type`` selects among the dataset types configured in the ``caom`` section, and reports an error if given a dataset type that has no CAOM configuration.

Installation
============

The command needs the ``caom2`` package, which is an optional dependency:

.. code-block:: bash

   pip install lsst-dax-obscore[caom]

Nothing else in ``dax_obscore`` imports ``caom2``, so an installation without it is fully functional for every other command.
Attempting to run ``export-caom`` without it raises an error naming the extra.

Configuration
=============

CAOM settings live in a ``caom`` section of the ordinary ObsCore configuration file.
Only dataset types listed under ``caom.dataset_types`` are exported, and each must also appear in the top-level ``dataset_types``.

.. code-block:: yaml

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
       visit_image:
         observation_id_fmt: "{records[visit].name}"
         product_id_fmt: "visit_image-{records[detector].full_name}"
         s_pixel_scale: 0.2
       deep_coadd:
         observation_id_fmt: "{skymap}-{tract}"
         product_id_fmt: "deep_coadd-{patch}-{band}"
         derived: true
         algorithm: "lsst.deep_coadd.DRP.DP1.DM-51335"
         s_pixel_scale: 0.2
         auxiliary_datasets:
           deep_coadd_n_image: auxiliary

The attributes are documented on :py:class:`~lsst.dax.obscore.caom_config.CaomConfig`, :py:class:`~lsst.dax.obscore.caom_config.CaomDatasetTypeConfig` and :py:class:`~lsst.dax.obscore.caom_config.CaomProvenanceConfig`.

Values derived from the ObsCore configuration
---------------------------------------------

Anything the ObsCore configuration already states is derived rather than repeated.

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - CAOM element
     - Source
   * - ``Observation.collection``
     - ``obs_collection``
   * - ``Observation.instrument``
     - record ``instrument_name``, which honors ``fallback_instrument``
   * - ``Observation.target``
     - record ``target_name``
   * - ``Observation.metaRelease``
     - ``origin.publication_date``
   * - ``Telescope.name``
     - ``caom.telescope_name``, defaulting to ``facility_name``
   * - ``Telescope.geoLocation{X,Y,Z}``
     - ``astropy.coordinates.EarthLocation.of_site`` applied to the facility name
   * - ``Proposal.title``
     - ``origin.title``
   * - ``Proposal.project``
     - ``obs_collection``
   * - ``Provenance.producer``
     - ``origin.publisher``
   * - ``Provenance.reference``
     - ``origin.reference_url``
   * - ``Provenance.runID``
     - the Butler RUN collection of the dataset
   * - ``Plane.calibrationLevel``
     - ``dataset_types.<X>.calib_level``
   * - ``Plane.dataProductType``
     - ``dataset_types.<X>.dataproduct_type``
   * - ``Plane.observable``
     - ``dataset_types.<X>.o_ucd``
   * - ``Position.bounds``
     - the spatial region held in the registry
   * - ``Position.dimension``
     - ``dataset_types.<X>.s_xel``
   * - ``Position.sampleSize``
     - ``caom.dataset_types.<X>.s_pixel_scale``, converted from arcsec to degrees
   * - ``Energy.bounds``
     - ``spectral_ranges``, through the record ``em_min`` and ``em_max``
   * - ``Energy.bandpassName``
     - record ``em_filter_name``
   * - ``Time.bounds`` and ``Time.exposure``
     - record ``t_min``, ``t_max`` and ``t_exptime``

``EarthLocation.of_site`` resolves ``Rubin:Simonyi``, so the ObsCore facility name is normally sufficient.
Astropy does not know every facility code, so ``caom.geo_location`` supplies the geocentric position directly when the lookup fails.

Notes on individual attributes
------------------------------

``provenance.name`` is the collection-specific common name of the process that produced the data, such as ``LSST Science Pipelines``.
It names the software, not the dataset, and CADC collections use values of the same kind.
``provenance.version`` is the version of that software, a fixed release-level value such as ``v29.0``.
The per-execution value is the Butler RUN collection, which is reported as ``Provenance.runID``.

``observation_type`` is the ``OBSTYPE``-style description of the exposure.
The ``exposure`` dimension carries an ``observation_type`` field with values such as ``science``, ``bias`` and ``dark``, so exposure-based dataset types should read it from the record with ``observation_type_fmt``.
The ``visit`` dimension has no such field, because a visit is always an on-sky science observation, so visit-based and coadd dataset types take the literal default.
That default is ``science`` rather than the CAOM-conventional ``OBJECT`` so that literal and templated values share one vocabulary.
``intent_fmt`` exists for the same reason: a ``raw`` bias is a calibration rather than a science observation.

``content_type`` varies by dataset type and is not always FITS.
It cannot reuse the ObsCore ``access_format``, which describes the DataLink response rather than the file.

``algorithm`` applies only when ``derived`` is true.
The ``caom2`` library forces a ``SimpleObservation`` to use the algorithm name ``exposure`` and requires a ``DerivedObservation`` to use anything else, so an ``algorithm`` on a non-derived dataset type is rejected, as is a derived dataset type without one.

Template attributes
===================

Every ``_fmt`` attribute is a :py:meth:`str.format` template expanded against one namespace holding:

- the names of the ``DataId`` attributes, for example ``{tract}`` or ``{detector}``,
- ``records`` indexed by dimension name, for example ``{records[visit].name}``,
- ``id``, the dataset UUID, ``run``, the RUN collection, and ``dataset_type``,
- ``butler_uri``, the resolved Butler URI of the file, which is an empty string when the repository has no datastore,
- every column of the ObsCore record, for example ``{obs_id}`` or ``{lsst_patch}``.

Auxiliary dataset types have no ObsCore record, so a template used for one may refer only to the data ID, ``id``, ``run``, ``dataset_type`` and ``butler_uri``.

Identifiers must be legal in a CAOM URI
---------------------------------------

CAOM identifies an observation by a ``caom:<collection>/<observation_id>`` URI, and a URI path component may not contain a space, slash, backslash or percent character.
An expanded ``observation_id_fmt`` or ``product_id_fmt`` containing one of those characters is reported as an error naming the template.
The identifier is never rewritten to make it legal, because it is the identifier CADC ingests; the configuration must produce a legal one.

A skymap named ``discrete/ci_hsc`` therefore cannot be used in an observation identifier, while ``lsst_cells_v1`` can.

Observation merging
===================

Dataset types whose ``observation_id_fmt`` templates expand to the same value contribute planes to a single Observation.
For a visit, ``raw`` contributes a calibration-level-1 plane, ``visit_image`` a level-2 plane and ``difference_image`` a level-3 plane, and CADC receives one ingestable document per visit.
This is the canonical CAOM shape, and it is the reason the exporter accumulates observations rather than streaming them.

The grouping is entirely a property of the templates.
For DP1 the shapes are one Observation per visit with one plane per detector for visit-based data, and one Observation per tract with one plane per patch and band for coadds.
A tract observation identifier and a visit observation identifier never collide, so coadds never merge with visit data.

Two rules apply when dataset types share an observation, and both are logged.

- The first dataset type to reach an observation identifier sets the observation-level attributes.
  A later dataset type that disagrees about the instrument, intent, type or ``derived`` flag produces a warning and is ignored, because there is no basis for choosing between them.
- Two dataset types producing the same ``(observation_id, product_id)`` pair abort the export.
  Dataset types sharing an observation are expected to carry distinct product identifiers.

Memory and large exports
========================

Every observation is held in memory until the export finishes, because a single observation can gather planes from several dataset types that are queried in separate passes.
At DP1 scale this is a few thousand observations and is not a concern.
Use ``--where`` to divide a larger export into chunks.

Known gaps
==========

This command is a prototype, and the following are recorded rather than solved.

- **The dataset DOI has no home in CAOM.**
  CAOM 2.4 has no DOI attribute.
  ``Plane.creatorID`` is an IVOA dataset identifier and is not a DOI, so ``origin.citation`` is not exported.
- **``Artifact.contentLength`` is not set.**
  The Butler file datastore records carry a ``file_size``, which needs no network access, but that value is sometimes ``-1``.
  The intended resolution is to read the datastore records and fall back to ``ResourcePath.size()`` for those cases behind an opt-in flag, since the fallback issues one HTTP ``HEAD`` per file.
- **There are no ``Part`` or ``Chunk`` elements.**
  CAOM 2.5 removes them, and the discovery metadata sits at plane level.
  A per-dataset WCS could be reconstructed from the registry region if they prove necessary before then, still without reading a file.
- **``s_pixel_scale`` is a placeholder.**
  It belongs beside ``s_xel`` in the ObsCore ``DatasetTypeConfig`` so that one statement of a dataset type's geometry serves both exports, and it moves there once ``daf_butler`` supports it.
  It is read through a single accessor, :py:meth:`~lsst.dax.obscore.caom_config.CaomConfig.pixel_scale_for`, so that the migration touches one function.
- **ObsCore ``s_resolution`` is never populated.**
  It is the PSF full width at half maximum, which is measured per dataset and is not available from registry records, so ``Position.resolution`` is left unset.
  It is deliberately not filled with the pixel scale, which is a different quantity.
