# coding=utf-8

# Copyright © 2016 Computational Molecular Biology Group,
#                  Freie Universität Berlin (GER)
#
# This file is part of ReaDDy.
#
# ReaDDy is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General
# Public License along with this program. If not, see
# <http://www.gnu.org/licenses/>.

import os
import unittest
import tempfile
import shutil
import h5py
import numpy as np
import readdy._internal.common as common
import readdy._internal.common.io as io


class TestSchemeApi(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.dir = tempfile.mkdtemp("test-io")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.dir, ignore_errors=True)

    def test_readwrite_double_and_string(self):
        fname = os.path.join(self.dir, "test_readwrite_double_and_string.h5")
        data = np.array([[2.222, 3, 4, 5], [3.3, 3, 3, 3]], dtype=np.float64)
        f = io.File(fname, io.FileAction.CREATE)
        f.write_double("/sowas", data)
        f.write_string("/maeh", u"hierstehtwas")
        f.close()

        with h5py.File(fname, "r") as f2:
            np.testing.assert_equal(f2.get('/sowas'), data)
            np.testing.assert_equal(f2.get("/maeh").value, u"hierstehtwas")

    def test_groups_readwrite(self):
        fname = os.path.join(self.dir, "test_groups_readwrite.h5")
        data = np.array([[2.222, 3, 4, 5], [3.3, 3, 3, 3]], dtype=np.float64)
        f = io.File(fname, io.FileAction.CREATE)
        g = f.create_group("/my_super_group")
        subg = g.create_group("my_super_subgroup")
        g.write_double("doubleds", data)
        subg.write_string("stringds", u"jap")
        f.close()

        with h5py.File(fname, "r") as f2:
            np.testing.assert_equal(f2.get("/my_super_group")["doubleds"], data)
            np.testing.assert_equal(f2.get("/my_super_group").get("my_super_subgroup")["stringds"].value, u"jap")
            np.testing.assert_equal(f2.get("/my_super_group/my_super_subgroup")["stringds"].value, u"jap")

    def test_append(self):
        fname = os.path.join(self.dir, "test_append.h5")
        f = io.File(fname, io.FileAction.CREATE)
        g = f.create_group("/append_group")

        ds_double = io.DataSet_double("doubleds", g, [io.DataSet_double.UNLIMITED_DIMS], [0])

        f.close()
        pass

if __name__ == '__main__':
    common.set_logging_level("debug")
    unittest.main()
