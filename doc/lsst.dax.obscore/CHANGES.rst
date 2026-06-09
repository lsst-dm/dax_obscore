lsst-dax-obscore v30.0.8 (2026-06-09)
=====================================

New Features
------------

- Added ``butler obscore export-regions`` subcommand that writes a FITS binary table of distinct IVOA STC-S region strings for datasets matching a dataset type, single collection, and optional where query.
  Suitable as input to the TOPCAT STILTS ``mocshape`` command. (`DM-53569 <https://rubinobs.atlassian.net/browse/DM-53569>`_)
- Added an SIAv2 configuration for Prompt Data Products. (`DM-54169 <https://rubinobs.atlassian.net/browse/DM-54169>`_)


lsst-dax-obscore v30.0.0 (2026-01-16)
=====================================

Changes to configurations
-------------------------

- Modified the DP1 ObsTAP configuration to match the currently expected dataset types. (`DM-49669 <https://rubinobs.atlassian.net/browse/DM-49669>`_)


Miscellaneous Changes of Minor Interest
---------------------------------------

- Added ``obs_publisher_did`` format string to configurations. (`DM-51383 <https://rubinobs.atlassian.net/browse/DM-51383>`_)
- * Removed detector from ``obs_id`` field.
  * Fixed data type from arrow schema for integer fields.
  * Corrected the data type for `int`` fields to "long" where needed in ObsCore model.
  * Added pixel dimensions to DP1 config. (`DM-51495 <https://rubinobs.atlassian.net/browse/DM-51495>`_)
- Added ``COOSYS`` and ``TIMESYS`` references to the VOTable fields. (`DM-51564 <https://rubinobs.atlassian.net/browse/DM-51564>`_)


lsst-dax-obscore v29.0.0 (2025-03-26)
=====================================

The package is now available on `PyPI <https://pypi.org/project/lsst-dax-obscore/>`_.

New Features
------------

* Tested with Python 3.13. (`DM-47128 <https://rubinobs.atlassian.net/brows/DM-47128>`_)
* Added ability for entry points to register overrides for the handling of SIAv2 queries to conform to alternate dimension universes. (`DM-48282 <https://rubinobs.atlassian.net/brows/DM-48282>`_)
* Modified the subcommand loading to use entry points.
  The ``DAF_BUTLER_PLUGINS`` environment variable is not longer set. (`DM-42226 <https://rubinobs.atlassian.net/brows/DM-42226>`_)

Performance Improvements
------------------------

* Reorganized the record exporter to be more efficient. (`DM-47980 <https://rubinobs.atlassian.net/brows/DM-47980>`_)

Configuration Changes
---------------------

* Added initial configuration for Data Preview 1. (`DM-49017 <https://rubinobs.atlassian.net/brows/DM-49017>`_)
* Added initial configuration for use of SIAv2 with embargo rack. (`DM-46990 <https://rubinobs.atlassian.net/brows/DM-46990>`_)
* Corrected ``access_format`` to refer to data link and not FITS. (`DM-47977 <https://rubinobs.atlassian.net/brows/DM-47977>`_)

lsst-dax-obscore v28.0.0 (2024-10-04)
=====================================

New Features
------------

* Added SIAv2-over-butler interface.
  It is now possible to parse SIAv2-style query parameters, query Butler, and return ObsCore records in VOTable format.
  Added ``butler obscore siav2`` subcommand to allow SIAv2-style queries to be run directly on a butler repository.
  The implementation was written up `for the ADASS conference <https://arxiv.org/abs/2501.00544>_`.
  (`DM-45860 <https://rubinobs.atlassian.net/brows/DM-45860>`_)

Configuration Changes
---------------------

* Improved DP0.2 configuration to more correctly reflect our standard data model. (`DM-34685 <https://rubinobs.atlassian.net/brows/DM-34685>`_)

Miscellaneous Changes of Minor Interest
---------------------------------------

* Changed the internal query system to support remote butler. (`DM-46363 <https://rubinobs.atlassian.net/brows/DM-46363>`_)


lsst-dax-obscore v27.0.0 (2024-01-08)
=====================================

First official release of the package tied to a science pipelines release.
Supports exporting of static ObsCore records from a butler repository. (`DM-34483 <https://rubinobs.atlassian.net/brows/DM-34483>`_)

Miscellaneous Changes of Minor Interest
---------------------------------------

* Reorganized the package layout and linting to match the current middleware standards and fix deprecation warnings. (`DM-40287 <https://rubinobs.atlassian.net/brows/DM-40287>`_)
