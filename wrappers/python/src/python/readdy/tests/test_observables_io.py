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


class TestObservablesIO(unittest.TestCase):
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
            handle.enable_write_to_file(f, u"particle_positions", int(3))
            sim.run_scheme_readdy(True).configure().run(n_timesteps)
            handle.flush()

        with h5py.File(fname, "r") as f2:
            data = f2["readdy/observables/particle_positions"][:]
            np.testing.assert_equal(len(data), n_timesteps + 1)
            for t, positions in enumerate(data):
                # we begin with two particles
                np.testing.assert_equal(len(positions), t + 2)
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
        typeid_A = sim.register_particle_type("A", .1, .1)
        typeid_B = sim.register_particle_type("B", .1, .1)
        sim.add_particle("A", common.Vec(0, 0, 0))
        sim.add_particle("B", common.Vec(0, 0, 0))
        # every time step, add one particle
        sim.register_observable_n_particles(1, lambda n: sim.add_particle("A", common.Vec(1.5, 2.5, 3.5)), ["A"])
        handle = sim.register_observable_particles(1, None)
        n_timesteps = 19
        with closing(io.File(fname, io.FileAction.CREATE, io.FileFlag.OVERWRITE)) as f:
            handle.enable_write_to_file(f, u"particles", int(3))
            sim.run_scheme_readdy(True).configure().run(n_timesteps)
            handle.flush()

        with h5py.File(fname, "r") as f2:
            types = f2["readdy/observables/particles/types"][:]
            ids = f2["readdy/observables/particles/ids"][:]
            positions = f2["readdy/observables/particles/positions"][:]
            for t in range(n_timesteps):
                np.testing.assert_equal(len(types[t]), t + 3)
                np.testing.assert_equal(len(ids[t]), t + 3)
                np.testing.assert_equal(len(positions[t]), t + 3)
                np.testing.assert_equal(types[t][0], typeid_A)
                np.testing.assert_equal(positions[t][0]["x"], 0)
                np.testing.assert_equal(positions[t][0]["y"], 0)
                np.testing.assert_equal(positions[t][0]["z"], 0)
                np.testing.assert_equal(positions[t][1]["x"], 0)
                np.testing.assert_equal(positions[t][1]["y"], 0)
                np.testing.assert_equal(positions[t][1]["z"], 0)
                np.testing.assert_equal(types[t][1], typeid_B)
                for others in range(2, len(types[t])):
                    np.testing.assert_equal(types[t][others], typeid_A)
                    np.testing.assert_equal(positions[t][others]["x"], 1.5)
                    np.testing.assert_equal(positions[t][others]["y"], 2.5)
                    np.testing.assert_equal(positions[t][others]["z"], 3.5)

    def test_radial_distribution_observable(self):
        common.set_logging_level("error")
        fname = os.path.join(self.dir, "test_observables_radial_distribution.h5")

        simulation = Simulation()
        simulation.set_kernel("SingleCPU")

        box_size = common.Vec(10, 10, 10)
        simulation.kbt = 2
        simulation.periodic_boundary = [True, True, True]
        simulation.box_size = box_size
        simulation.register_particle_type("A", .2, 1.)
        simulation.register_particle_type("B", .2, 1.)
        simulation.register_potential_harmonic_repulsion("A", "B", 10)
        simulation.add_particle("A", common.Vec(-2.5, 0, 0))
        simulation.add_particle("B", common.Vec(0, 0, 0))
        bin_borders = np.arange(0, 5, .01)
        density = 1. / (box_size[0] * box_size[1] * box_size[2])
        n_time_steps = 50
        callback_centers = []
        callback_rdf = []

        def rdf_callback(pair):
            callback_centers.append(pair[0])
            callback_rdf.append(pair[1])

        handle = simulation.register_observable_radial_distribution(1, rdf_callback, bin_borders, ["A"], ["B"], density)
        with closing(io.File(fname, io.FileAction.CREATE, io.FileFlag.OVERWRITE)) as f:
            handle.enable_write_to_file(f, u"radial_distribution", int(3))
            simulation.run(n_time_steps, 0.02)
            handle.flush()

        with h5py.File(fname, "r") as f2:
            bin_centers = f2["readdy/observables/radial_distribution/bin_centers"][:]
            distribution = f2["readdy/observables/radial_distribution/distribution"][:]
            for t in range(n_time_steps):
                np.testing.assert_equal(bin_centers, np.array(callback_centers[t]))
                np.testing.assert_equal(distribution[t], np.array(callback_rdf[t]))

    def test_center_of_mass_observable(self):
        common.set_logging_level("error")
        fname = os.path.join(self.dir, "test_observables_com.h5")

        simulation = Simulation()
        simulation.set_kernel("SingleCPU")

        box_size = common.Vec(10, 10, 10)
        simulation.kbt = 2
        simulation.periodic_boundary = [True, True, True]
        simulation.box_size = box_size
        simulation.register_particle_type("A", .2, 1.)
        simulation.register_particle_type("B", .2, 1.)
        simulation.add_particle("A", common.Vec(-2.5, 0, 0))
        simulation.add_particle("B", common.Vec(0, 0, 0))
        n_time_steps = 50
        callback_com = []

        def com_callback(vec):
            callback_com.append(vec)

        handle = simulation.register_observable_center_of_mass(1, com_callback, ["A", "B"])
        with closing(io.File(fname, io.FileAction.CREATE, io.FileFlag.OVERWRITE)) as f:
            handle.enable_write_to_file(f, u"com", 3)
            simulation.run(n_time_steps, 0.02)
            handle.flush()

        with h5py.File(fname, "r") as f2:
            com = f2["readdy/observables/com"][:]
            for t in range(n_time_steps):
                np.testing.assert_equal(com[t]["x"], callback_com[t][0])
                np.testing.assert_equal(com[t]["y"], callback_com[t][1])
                np.testing.assert_equal(com[t]["z"], callback_com[t][2])

    def test_histogram_along_axis_observable(self):
        common.set_logging_level("error")
        fname = os.path.join(self.dir, "test_observables_hist_along_axis.h5")

        simulation = Simulation()
        simulation.set_kernel("SingleCPU")

        box_size = common.Vec(10, 10, 10)
        simulation.kbt = 2
        simulation.periodic_boundary = [True, True, True]
        simulation.box_size = box_size
        simulation.register_particle_type("A", .2, 1.)
        simulation.register_particle_type("B", .2, 1.)
        simulation.register_potential_harmonic_repulsion("A", "B", 10)
        simulation.add_particle("A", common.Vec(-2.5, 0, 0))
        simulation.add_particle("B", common.Vec(0, 0, 0))
        bin_borders = np.arange(0, 5, .01)
        n_time_steps = 50
        callback_hist = []

        def hist_callback(hist):
            callback_hist.append(hist)

        handle = simulation.register_observable_histogram_along_axis(1, hist_callback, bin_borders, 0, ["A", "B"])
        with closing(io.File(fname, io.FileAction.CREATE, io.FileFlag.OVERWRITE)) as f:
            handle.enable_write_to_file(f, u"hist_along_x_axis", int(3))
            simulation.run(n_time_steps, 0.02)
            handle.flush()

        with h5py.File(fname, "r") as f2:
            histogram = f2["readdy/observables/hist_along_x_axis"][:]
            for t in range(n_time_steps):
                np.testing.assert_equal(histogram[t], np.array(callback_hist[t]))

    def test_write_multiple_observables_at_once(self):
        # todo
        pass

if __name__ == '__main__':
    unittest.main()
