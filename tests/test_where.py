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

from lsst.dax.obscore.config import WhereBind


class WhereBindTestCase(unittest.TestCase):
    """Test the WhereBind model."""

    def test_where(self):
        w1 = WhereBind(where="instrument = 'LSSTCam'")
        w2 = WhereBind(where="physical_filter IN (phys)", bind={"phys": ("filter1", "filter2")})

        w_and = WhereBind.combine([w1, w2])
        w_ref = WhereBind(
            where="(instrument = 'LSSTCam') AND (physical_filter IN (phys))",
            bind={"phys": ("filter1", "filter2")},
        )
        self.assertEqual(w_and, w_ref)

        w_or = WhereBind.combine([w1, w2], mode="OR")
        w_ref = WhereBind(
            where="(instrument = 'LSSTCam') OR (physical_filter IN (phys))",
            bind={"phys": ("filter1", "filter2")},
        )
        self.assertEqual(w_or, w_ref)

        # Check that parentheses look okay combining OR and AND
        w_and_or = WhereBind.combine([w_or, w_and])
        w_ref = WhereBind(
            where="((instrument = 'LSSTCam') OR (physical_filter IN (phys))) "
            "AND ((instrument = 'LSSTCam') AND (physical_filter IN (phys)))",
            bind={"phys": ("filter1", "filter2")},
        )
        self.assertEqual(w_and_or, w_ref)

        # Duplicate bind keys but identical.
        w_dupe = WhereBind.combine([w2, w2])
        self.assertEqual(
            w_dupe,
            WhereBind(
                where="(physical_filter IN (phys)) AND (physical_filter IN (phys))",
                bind={"phys": ("filter1", "filter2")},
            ),
        )

        # Duplicates are ignored if there is no where string at all.
        w_empty = WhereBind(bind={"phys": ("filter1", "filter4")}, extras={"exposure"})
        w_and_empty = WhereBind.combine([w2, w_empty])
        self.assertEqual(
            w_and_empty,
            WhereBind(
                where="(physical_filter IN (phys))",
                bind={"phys": ("filter1", "filter2")},
                extras={"exposure"},
            ),
        )

        # Reuse phys bind key with different value and add instrument key.
        w3 = WhereBind(
            where="physical_filter IN (phys) AND instrument = INST",
            bind={"phys": ("filter0", "filter2"), "INST": "LATISS"},
        )

        with self.assertRaises(ValueError):
            WhereBind.combine([w2, w3])

        # Same as w3 but one duplicate key is identical and one is different.
        w4 = WhereBind(
            where="physical_filter IN (phys) AND instrument = INST",
            bind={"phys": ("filter0", "filter2"), "INST": "LSSTCam"},
        )

        with self.assertRaises(ValueError):
            WhereBind.combine([w1, w3, w4])

        with self.assertRaises(ValueError):
            WhereBind.combine([])


if __name__ == "__main__":
    unittest.main()
