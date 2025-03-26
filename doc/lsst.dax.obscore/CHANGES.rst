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
