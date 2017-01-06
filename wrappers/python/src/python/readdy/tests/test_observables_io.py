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
import readdy._internal.common as common
import readdy._internal.common.io as io
from readdy._internal.api import Simulation
from contextlib import closing
import numpy as np


class TestSchemeApi(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.dir = tempfile.mkdtemp("test-observables-io")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.dir, ignore_errors=True)

    def test_particle_positions_observable(self):
        fname = os.path.join(self.dir, "test_observables_particle_positions.h5")
        sim = Simulation()
        sim.set_kernel("SingleCPU")
        sim.box_size = common.Vec(13, 13, 13)
        sim.register_particle_type("A", .1, .1)
        sim.add_particle("A", common.Vec(0, 0, 0))
        # every time step, add one particle
        sim.register_observable_n_particles(1, lambda n: sim.add_particle("A", common.Vec(1.5, 2.5, 3.5)), ["A"])
        handle = sim.register_observable_particle_positions(1, None, [])
        n_timesteps = 19
        with closing(io.File(fname, io.FileAction.CREATE, io.FileFlag.OVERWRITE)) as f:
            handle.enable_write_to_file(f, "particle_positions", 3)
            sim.run_scheme_readdy(True).configure().run(n_timesteps)
            handle.flush()

        with h5py.File(fname, "r") as f2:
            data = f2["readdy/observables/particle_positions"][:]
            np.testing.assert_equal(len(data), n_timesteps+1)
            for t, positions in enumerate(data):
                # we begin with two particles
                np.testing.assert_equal(len(positions), t+2)
                np.testing.assert_equal(positions[0]["x"], 0)
                np.testing.assert_equal(positions[0]["y"], 0)
                np.testing.assert_equal(positions[0]["z"], 0)
                for i in range(1, len(positions)):
                    np.testing.assert_equal(positions[i]["x"], 1.5)
                    np.testing.assert_equal(positions[i]["y"], 2.5)
                    np.testing.assert_equal(positions[i]["z"], 3.5)
        common.set_logging_level("error")

    def test_particles_observable(self):
        fname = os.path.join(self.dir, "test_observables_particles.h5")
        sim = Simulation()
        sim.set_kernel("SingleCPU")
        sim.box_size = common.Vec(13, 13, 13)
        sim.register_particle_type("A", .1, .1)
        sim.register_particle_type("B", .1, .1)
        sim.add_particle("A", common.Vec(0, 0, 0))
        sim.add_particle("B", common.Vec(0, 0, 0))
        # every time step, add one particle
        sim.register_observable_(1, lambda n: sim.add_particle("A", common.Vec(1.5, 2.5, 3.5)), ["A"])
        handle = sim.register_observable_particle_positions(1, None, [])
        n_timesteps = 19
        with closing(io.File(fname, io.FileAction.CREATE, io.FileFlag.OVERWRITE)) as f:
            handle.enable_write_to_file(f, "particle_positions", 3)
            sim.run_scheme_readdy(True).configure().run(n_timesteps)
            handle.flush()

        with h5py.File(fname, "r") as f2:
            data = f2["readdy/observables/particle_positions"][:]
            np.testing.assert_equal(len(data), n_timesteps+1)
            for t, positions in enumerate(data):
                # we begin with two particles
                np.testing.assert_equal(len(positions), t+2)
                np.testing.assert_equal(positions[0]["x"], 0)
                np.testing.assert_equal(positions[0]["y"], 0)
                np.testing.assert_equal(positions[0]["z"], 0)
                for i in range(1, len(positions)):
                    np.testing.assert_equal(positions[i]["x"], 1.5)
                    np.testing.assert_equal(positions[i]["y"], 2.5)
                    np.testing.assert_equal(positions[i]["z"], 3.5)
        common.set_logging_level("error")


if __name__ == '__main__':
    unittest.main()
