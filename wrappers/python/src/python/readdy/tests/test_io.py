# coding=utf-8

# Copyright © 2018 Computational Molecular Biology Group,
#                  Freie Universität Berlin (GER)
#
# Redistribution and use in source and binary forms, with or
# without modification, are permitted provided that the
# following conditions are met:
#  1. Redistributions of source code must retain the above
#     copyright notice, this list of conditions and the
#     following disclaimer.
#  2. Redistributions in binary form must reproduce the above
#     copyright notice, this list of conditions and the following
#     disclaimer in the documentation and/or other materials
#     provided with the distribution.
#  3. Neither the name of the copyright holder nor the names of
#     its contributors may be used to endorse or promote products
#     derived from this software without specific
#     prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
# CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
# INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
# STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
# ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""
@author: clonker
"""

import os
import unittest
import tempfile
import shutil
import h5py
import numpy as np
import readdy._internal.readdybinding.common as common
import readdy._internal.readdybinding.common.io as io
from readdy._internal.readdybinding.api import Simulation

from readdy.util.testing_utils import ReaDDyTestCase
from readdy.util.trajectory_utils import TrajectoryReader
from contextlib import closing


class TestSchemeApi(ReaDDyTestCase):
    @classmethod
    def setUpClass(cls):
        cls.dir = tempfile.mkdtemp("test-io")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.dir, ignore_errors=True)

    def test_write_trajectory(self):
        traj_fname = os.path.join(self.dir, "traj.h5")
        simulation = Simulation("SingleCPU")
        simulation.context.box_size = [5., 5., 5.]
        simulation.context.particle_types.add("A", 0.0)
        simulation.context.reactions.add_conversion("A->A", "A", "A", 1.)

        def callback(_):
            simulation.add_particle("A", common.Vec(0, 0, 0))

        simulation.register_observable_n_particles(1, ["A"], callback)
        traj_handle = simulation.register_observable_trajectory(0)
        with closing(io.File.create(traj_fname, io.FileFlag.OVERWRITE)) as f:
            traj_handle.enable_write_to_file(f, u"", 3)
            loop = simulation.create_loop(1.)
            loop.write_config_to_file(f)
            loop.run(20)

        r = TrajectoryReader(traj_fname)
        trajectory_items = r[:]
        for idx, items in enumerate(trajectory_items):
            np.testing.assert_equal(len(items), idx+1)
            for item in items:
                np.testing.assert_equal(item.t, idx)
                np.testing.assert_equal(item.position, np.array([.0, .0, .0]))

        with h5py.File(traj_fname) as f:
            np.testing.assert_equal("A", f["readdy/config/particle_types"][0]["name"])

    def test_write_flat_trajectory(self):
        import readdy
        rds = readdy.ReactionDiffusionSystem([5, 5, 5])
        rds.add_species("A", 0.)
        simulation = rds.simulation()
        simulation.output_file = os.path.join(self.dir, "flat_traj.h5")

        def callback(_):
            simulation.add_particle("A", common.Vec(0, 0, 0))
        simulation.observe.number_of_particles(1, ["A"], callback=callback)
        simulation.record_trajectory(1)
        simulation.run(20, 1, show_system=False)

        r = TrajectoryReader(simulation.output_file)
        trajectory_items = r[:]
        for idx, items in enumerate(trajectory_items):
            np.testing.assert_equal(len(items), idx+1)
            for item in items:
                np.testing.assert_equal(item.t, idx)
                np.testing.assert_equal(item.position, np.array([.0, .0, .0]))

    def test_write_trajectory_as_observable(self):
        traj_fname = os.path.join(self.dir, "traj_as_obs.h5")
        simulation = Simulation("SingleCPU")
        simulation.context.box_size = [5., 5., 5.]
        simulation.context.particle_types.add("A", 0.0)

        def callback(_):
            simulation.add_particle("A", common.Vec(0, 0, 0))

        simulation.register_observable_n_particles(1, ["A"], callback)
        traj_handle = simulation.register_observable_trajectory(1)

        with closing(io.File.create(traj_fname, io.FileFlag.OVERWRITE)) as f:
            traj_handle.enable_write_to_file(f, u"", int(3))
            simulation.run(20, 1)

        r = TrajectoryReader(traj_fname)
        trajectory_items = r[:]
        for idx, items in enumerate(trajectory_items):
            np.testing.assert_equal(len(items), idx+1)
            for item in items:
                np.testing.assert_equal(item.t, idx)
                np.testing.assert_equal(item.position, np.array([.0, .0, .0]))


    def test_open_and_discover_file(self):
        fname = os.path.join(self.dir, "test_open_and_discover_file.h5")
        data = np.array([[2.222, 3, 4, 5], [3.3, 3, 3, 3]], dtype=np.float64)
        with closing(io.File.create(fname)) as f:
            g = f.create_group("/my_super_group")
            subg = g.create_group("my_super_subgroup")
            g.write_double("doubleds", data)
            subg.write_string("stringds", u"jap")

        with closing(io.File.open(fname, flag=io.FileFlag.READ_WRITE)) as f:
            sg = f.subgroups()
            np.testing.assert_equal(len(sg), 1)
            np.testing.assert_equal(sg[0], u"my_super_group")
            sub_group = f.get_subgroup(sg[0])
            np.testing.assert_equal(len(sub_group.data_sets()), 1)
            np.testing.assert_equal(sub_group.data_sets()[0], "doubleds")
            subsub_groups = sub_group.subgroups()
            np.testing.assert_equal(len(subsub_groups), 1)
            np.testing.assert_equal(subsub_groups[0], "my_super_subgroup")
            super_subgroup = sub_group.get_subgroup("my_super_subgroup")
            np.testing.assert_equal(len(super_subgroup.data_sets()), 1)
            np.testing.assert_equal(super_subgroup.data_sets()[0], "stringds")

    def test_readwrite_double_and_string(self):
        fname = os.path.join(self.dir, "test_readwrite_double_and_string.h5")
        data = np.array([[2.222, 3, 4, 5], [3.3, 3, 3, 3]], dtype=np.float64)
        with closing(io.File.create(fname)) as f:
            f.write_double("/sowas", data)
            f.write_string("/maeh", u"hierstehtwas")

        with h5py.File(fname, "r") as f2:
            np.testing.assert_equal(f2.get('/sowas'), data)
            np.testing.assert_equal(f2.get("/maeh").value.decode(), u"hierstehtwas")

    def test_groups_readwrite(self):
        fname = os.path.join(self.dir, "test_groups_readwrite.h5")
        data = np.array([[2.222, 3, 4, 5], [3.3, 3, 3, 3]], dtype=np.float64)
        with closing(io.File.create(fname)) as f:
            g = f.create_group("/my_super_group")
            subg = g.create_group("my_super_subgroup")
            g.write_double("doubleds", data)
            subg.write_string("stringds", u"jap")

        with h5py.File(fname, "r") as f2:
            np.testing.assert_equal(f2.get("/my_super_group")["doubleds"], data)
            np.testing.assert_equal(f2.get("/my_super_group").get("my_super_subgroup")["stringds"].value.decode(), u"jap")
            np.testing.assert_equal(f2.get("/my_super_group/my_super_subgroup")["stringds"].value.decode(), u"jap")


if __name__ == '__main__':
    unittest.main()
