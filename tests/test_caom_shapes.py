# This file is part of dax_obscore.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (http://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

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
