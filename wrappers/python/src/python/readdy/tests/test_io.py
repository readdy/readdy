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

"""
@author: clonker
"""

import os
import unittest
import tempfile
import shutil
import h5py
import numpy as np
import readdy._internal.common as common
import readdy._internal.common.io as io
from readdy._internal.api import Simulation
from readdy.util.trajectory_utils import TrajectoryReader


class TestSchemeApi(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.dir = tempfile.mkdtemp("test-io")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.dir, ignore_errors=True)

    def test_write_trajectory(self):
        common.set_logging_level("debug")
        traj_fname = os.path.join(self.dir, "traj.h5")
        simulation = Simulation()
        simulation.set_kernel("SingleCPU")
        simulation.set_time_step(1)
        simulation.box_size = common.Vec(5,5,5)
        simulation.register_particle_type("A", 0.0, 0.0)

        def callback(_):
            simulation.add_particle("A", common.Vec(0, 0, 0))

        simulation.register_observable_n_particles(1, callback, ["A"])
        simulation.record_trajectory(traj_fname, 0, 3)
        simulation.run_scheme_readdy(True).configure().run(20)
        simulation.close_trajectory_file()

        r = TrajectoryReader(traj_fname)
        trajectory_items = r[:]
        for idx, items in enumerate(trajectory_items):
            np.testing.assert_equal(len(items), idx+1)
            for item in items:
                np.testing.assert_equal(item.t, idx)
                np.testing.assert_equal(item.position, np.array([.0, .0, .0]))

        common.set_logging_level("debug")

    def test_readwrite_double_and_string(self):
        fname = os.path.join(self.dir, "test_readwrite_double_and_string.h5")
        data = np.array([[2.222, 3, 4, 5], [3.3, 3, 3, 3]], dtype=np.float64)
        f = io.File(fname, io.FileAction.CREATE)
        f.write_double("/sowas", data)
        f.write_string("/maeh", u"hierstehtwas")
        f.close()

        with h5py.File(fname, "r") as f2:
            np.testing.assert_equal(f2.get('/sowas'), data)
            np.testing.assert_equal(f2.get("/maeh").value.decode(), u"hierstehtwas")

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
            np.testing.assert_equal(f2.get("/my_super_group").get("my_super_subgroup")["stringds"].value.decode(), u"jap")
            np.testing.assert_equal(f2.get("/my_super_group/my_super_subgroup")["stringds"].value.decode(), u"jap")

    def test_append(self):
        fname = os.path.join(self.dir, "test_append.h5")
        f = io.File(fname, io.FileAction.CREATE)
        g = f.create_group("/append_group")

        full_data = [[3.3, 2.2], [1, 1], [3.4, 2.4], [14, 14], [5.5, 5.5]]

        ds_double = io.DataSet_double("doubleds", g, [2, 2], [io.unlimited_dims(), 2])
        ds_double.append(np.array([[3.3, 2.2], [1, 1]], dtype=np.float64))
        ds_double.append(np.array([[3.4, 2.4], [14, 14], [5.5, 5.5]], dtype=np.float64))

        f.close()

        with h5py.File(fname, "r") as f2:
            np.testing.assert_equal(f2.get("append_group")["doubleds"][:], full_data)


if __name__ == '__main__':
    common.set_logging_level("debug")
    unittest.main()
