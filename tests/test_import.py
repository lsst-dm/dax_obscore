# This file is part of dax_obscore.
#
# LSST Data Management System
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
# See COPYRIGHT file at the top of the source tree.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <https://www.lsstcorp.org/LegalNotices/>.

import unittest

from lsst.utils.tests import ImportTestCase


class PipeTasksImportTestCase(ImportTestCase):
    """Test that every file can be imported.

    dax_obscore does not import the command-line code at this time.
    """

    PACKAGES = {
        "lsst.dax.obscore.cli.cmd",
        "lsst.dax.obscore.script",
    }


if __name__ == "__main__":
    unittest.main()
