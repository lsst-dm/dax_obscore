.. py:currentmodule:: lsst.dax.obscore

.. _lsst.dax.obscore:

################
lsst.dax.obscore
################

This package implements ``butler`` CLI plugin ``obscore-export`` which generates ``ObsCore`` data model from Butler registry.

.. _lsst.dax.obscore-cli:

.. click:: lsst.daf.butler.cli.butler:cli
   :prog: butler
   :nested: full
   :commands: obscore-export


.. _lsst.dax.obscore-config:

Configuration
=============

``butler obscore-export`` requires a configuration file in YAML format, the data in that file is used to populate in-memory configuration consisting of instances of :py:class:`~lsst.dax.obscore.ExporterConfig` and :py:class:`~lsst.dax.obscore.DatasetTypeConfig` classes.
Documentation for these classes explains the type and meaning of the attributes.
An example configuration is defined in ``configs/example.yaml`` file.

Some attributes, such as ``obs_if_fmt``, represent format strings used to generate attribute values based on values of other items.
The syntax of these format strings corresponds to the syntax used by :py:meth:`str.format` method.
Possible list of attributes that can be used in the format strings includes:

- names of the ``DataId`` attributes, e.g. ``instrument``,
- ``records`` with a corresponding dimension name as an index and and attribute of that record. e.g. ``records[exposure].obs_id``,
- names of other attributes of ObsCore data model, e.g. ``dataproduct_type``, the list of attributes is limited to a subset, currently implemented in this module.

Few examples of specifying ``obs_id_fmt`` value::

   obs_id_fmt: "{skymap}-{tract}-{patch}"
   obs_id_fmt: "{records[exposure].obs_id}"
   obs_id_fmt: "{records[visit].name}"

``obscore-export`` will read datasets from Registry using collections specified in ``collections`` configuration attribute and dataset types that appear in ``dataset_types`` attribute (indexed by dataset type names).



.. _lsst.dax.obscore-pyapi:

Python API reference
====================

.. py:module:: lsst.dax.obscore

.. autoclass:: lsst.dax.obscore.DatasetTypeConfig
   :members:

.. autoclass:: lsst.dax.obscore.ExporterConfig
   :members:

.. autoclass:: lsst.dax.obscore.ObscoreExporter
   :members:
