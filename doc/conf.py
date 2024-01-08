"""Sphinx configuration file for an LSST stack package.

This configuration only affects single-package Sphinx documentation builds.
"""

from documenteer.conf.pipelinespkg import *  # noqa: F403, import *

project = "dax_obscore"
html_theme_options["logotext"] = project  # noqa: F405, unknown name
html_title = project
html_short_title = project
doxylink = {}
exclude_patterns = ["changes/*"]

intersphinx_mapping["pydantic"] = ("https://docs.pydantic.dev/latest/", None)  # noqa: F405
intersphinx_mapping["lsst"] = ("https://pipelines.lsst.io/v/weekly/", None)  # noqa: F405
