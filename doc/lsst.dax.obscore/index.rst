.. py:currentmodule:: lsst.dax.obscore

.. _lsst.dax.obscore:

################
lsst.dax.obscore
################

This package implements ``butler`` CLI plugin ``obscore-export`` which generates ``ObsCore`` data model from Butler registry.

.. _lsst.dax.obscore-changes:

Changes
=======

.. toctree::
   :maxdepth: 1

   CHANGES.rst

.. _lsst.dax.obscore-cli:

.. click:: lsst.daf.butler.cli.butler:cli
   :prog: butler
   :nested: full
   :commands: obscore-export


.. _lsst.dax.obscore-config:

Configuration
=============

``butler obscore-export`` requires a configuration file in YAML format, the data in that file is used to populate in-memory configuration consisting of instances of :py:class:`~lsst.dax.obscore.ExporterConfig` and :py:class:`~lsst.daf.butler.registry.obscore.DatasetTypeConfig` classes.
Documentation for these classes explains the type and meaning of the attributes.
An example configuration is defined in ``configs/example.yaml`` file.

Some attributes, such as ``obs_id_fmt``, represent format strings used to generate attribute values based on values of other items.
The syntax of these format strings corresponds to the syntax used by :py:meth:`str.format` method.
Possible list of attributes that can be used in the format strings includes:

- names of the ``DataId`` attributes, e.g. ``instrument``,
- ``records`` with a corresponding dimension name as an index and an attribute of that record. e.g. ``records[exposure].obs_id``,
- ``id`` which represents dataset ID (UUID),
- names of other attributes of ObsCore data model, e.g. ``dataproduct_type``, the list of attributes is limited to a subset currently implemented in this module.

Few examples of specifying ``obs_id_fmt`` value::

   obs_id_fmt: "{skymap}-{tract}-{patch}"
   obs_id_fmt: "{records[exposure].obs_id}"
   obs_id_fmt: "{records[visit].name}"

``extra_columns`` attributes can be used to provide additional global or per-dataset type columns with fixed values.
``extra_columns`` can also be used to override values for existing columns, e.g. override for an ``instrument_name`` can be specified as::

   extra_columns:
     instrument_name: LSSTCam

If global and per-dataset extra_column contain the same key, per-dataset value takes priority.

It is possible to use format strings in ``extra_columns`` values by providing ``template`` and ``type`` attributes instead of plain value, e.g.::

   extra_columns:
     day_obs:
       template: "{records[exposure].day_obs}"
       type: int

Supported data types are "int", "float", "str", and "bool".
If data type is not specified then "str" is assumed.
The "bool" type expects values (after template expansion) to be integer numbers, "0" maps to false, other values map to true.
If the attribute in the template string does not exist for a particular record then the whole value will be replaced with ``None/NULL``.

``obscore-export`` will read datasets from Registry using collections specified in ``collections`` configuration attribute and dataset types that appear in ``dataset_types`` attribute (indexed by dataset type names).


.. _lsst.dax.obscore-entry_points:

Entry Points
============

By default the SIAv2 query handler works with the default Butler dimension universe that is named ``daf_butler``.
Other dimension universes can be supported by providing a subclass to ``lsst.dax.obscore.siav2.SIAv2Handler`` specified in an entry point group named ``dax_obscore.siav2``.
The entry point label should be the dimension universe namespace.

For example, the default namespace entry point is defined as:

.. code:: toml

    [project.entry-points.'dax_obscore.siav2']
    daf_butler = "lsst.dax.obscore.siav2:get_daf_butler_siav2_handler"


.. _lsst.dax.obscore-pyapi:

Python API reference
====================

.. py:module:: lsst.dax.obscore

.. autoclass:: lsst.dax.obscore.ExporterConfig
   :members:

.. autoclass:: lsst.dax.obscore.ObscoreExporter
   :members:
