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
@author: chrisfroe
"""

from __future__ import print_function

import os
import unittest
import tempfile
import shutil
import h5py
import readdy._internal.readdybinding.common as common
import readdy._internal.readdybinding.common.io as io
from readdy._internal.readdybinding.api import Simulation
from readdy._internal.readdybinding.api import Context
from readdy._internal.readdybinding.api import KernelProvider

import readdy
from readdy.util import platform_utils
import readdy.util.io_utils as ioutils

from contextlib import closing
import numpy as np

from readdy.util.testing_utils import ReaDDyTestCase


class TestObservablesIO(ReaDDyTestCase):
    @classmethod
    def setUpClass(cls):
        cls.kernel_provider = KernelProvider.get()
        cls.kernel_provider.load_from_dir(str(platform_utils.get_readdy_plugin_dir()))
        cls.dir = tempfile.mkdtemp("test-observables-io")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.dir, ignore_errors=True)

    def test_particle_positions_observable(self):
        fname = os.path.join(self.dir, "test_observables_particle_positions.h5")
        context = Context()
        context.box_size = [13., 13., 13.]
        context.particle_types.add("A", .1)
        sim = Simulation("SingleCPU", context)
        sim.add_particle("A", common.Vec(0, 0, 0))
        # every time step, add one particle
        sim.register_observable_n_particles(1, ["A"], lambda n: sim.add_particle("A", common.Vec(1.5, 2.5, 3.5)))
        handle = sim.register_observable_particle_positions(1, [])
        n_timesteps = 19
        with closing(io.File.create(fname)) as f:
            handle.enable_write_to_file(f, u"particle_positions", int(3))
            sim.run(n_timesteps, 0)
            handle.flush()

        with h5py.File(fname, "r") as f2:
            data = f2["readdy/observables/particle_positions/data"][:]
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

    def test_virial_observable_CPU(self):
        fname = os.path.join(self.dir, "test_observables_virial.h5")
        context = Context()
        context.box_size = [13., 13., 13.]
        context.particle_types.add("A", .1)
        context.potentials.add_harmonic_repulsion("A", "A", 10., .5)
        sim = Simulation("CPU", context)
        for _ in range(10000):
            pos = common.Vec(*(13*np.random.random(size=3)-.5*13))
            sim.add_particle("A", pos)

        virials = []
        def virial_callback(virial):
            virials.append(np.ndarray((3,3), buffer=virial))

        handle = sim.register_observable_virial(1, virial_callback)
        with closing(io.File.create(fname)) as f:
            handle.enable_write_to_file(f, u"virial", int(3))
            sim.run(10, .1)
            handle.flush()

        with h5py.File(fname, "r") as f2:
            h5virials = f2["readdy/observables/virial/data"]
            for v, v2 in zip(virials, h5virials):
                np.testing.assert_almost_equal(v, v2)

    def test_virial_observable_SCPU(self):
        fname = os.path.join(self.dir, "test_observables_virial_scpu.h5")
        context = Context()
        context.box_size = [13., 13., 13.]
        context.particle_types.add("A", .1)
        context.potentials.add_harmonic_repulsion("A", "A", 10., .5)
        sim = Simulation("SingleCPU", context)
        for _ in range(10000):
            pos = common.Vec(*(13*np.random.random(size=3)-.5*13))
            sim.add_particle("A", pos)

        virials = []
        def virial_callback(virial):
            virials.append(np.ndarray((3,3), buffer=virial))

        handle = sim.register_observable_virial(1, virial_callback)
        with closing(io.File.create(fname)) as f:
            handle.enable_write_to_file(f, u"virial", int(3))
            sim.run(10, .1)
            handle.flush()

        with h5py.File(fname, "r") as f2:
            h5virials = f2["readdy/observables/virial/data"]
            for v, v2 in zip(virials, h5virials):
                np.testing.assert_almost_equal(v, v2)

    def test_particles_observable(self):
        fname = os.path.join(self.dir, "test_observables_particles.h5")
        context = Context()
        context.box_size = [13.,13.,13.]
        context.particle_types.add("A", .1)
        context.particle_types.add("B", .1)
        sim = Simulation("SingleCPU", context)
        sim.add_particle("A", common.Vec(0, 0, 0))
        sim.add_particle("B", common.Vec(0, 0, 0))
        # every time step, add one particle
        sim.register_observable_n_particles(1, ["A"], lambda n: sim.add_particle("A", common.Vec(1.5, 2.5, 3.5)))
        handle = sim.register_observable_particles(1)
        n_timesteps = 19
        with closing(io.File.create(fname)) as f:
            handle.enable_write_to_file(f, u"particles", int(3))
            loop = sim.create_loop(0)
            loop.write_config_to_file(f)
            loop.run(n_timesteps)
            handle.flush()

        from readdy.util.io_utils import get_particle_types

        particle_types = get_particle_types(fname)

        with h5py.File(fname, "r") as f2:

            types = f2["readdy/observables/particles/types"][:]
            ids = f2["readdy/observables/particles/ids"][:]
            positions = f2["readdy/observables/particles/positions"][:]
            for t in range(n_timesteps):
                np.testing.assert_equal(len(types[t]), t + 3)
                np.testing.assert_equal(len(ids[t]), t + 3)
                np.testing.assert_equal(len(positions[t]), t + 3)
                np.testing.assert_equal(types[t][0], particle_types["A"]["type_id"])
                np.testing.assert_equal(positions[t][0][0], 0)
                np.testing.assert_equal(positions[t][0][1], 0)
                np.testing.assert_equal(positions[t][0][2], 0)
                np.testing.assert_equal(positions[t][1][0], 0)
                np.testing.assert_equal(positions[t][1][1], 0)
                np.testing.assert_equal(positions[t][1][2], 0)
                np.testing.assert_equal(types[t][1], particle_types["B"]["type_id"])
                for others in range(2, len(types[t])):
                    np.testing.assert_equal(types[t][others], particle_types["A"]["type_id"])
                    np.testing.assert_equal(positions[t][others][0], 1.5)
                    np.testing.assert_equal(positions[t][others][1], 2.5)
                    np.testing.assert_equal(positions[t][others][2], 3.5)

    def test_radial_distribution_observable(self):
        fname = os.path.join(self.dir, "test_observables_radial_distribution.h5")

        context = Context()
        context.kbt = 2
        context.pbc = [True, True, True]
        box_size = [10.,10.,10.]
        context.box_size = box_size
        context.particle_types.add("A", .2)
        context.particle_types.add("B", .2)
        context.potentials.add_harmonic_repulsion("A", "B", 10, 2.)

        simulation = Simulation("SingleCPU", context)

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

        handle = simulation.register_observable_radial_distribution(1, bin_borders, ["A"], ["B"], density, rdf_callback)
        with closing(io.File.create(fname)) as f:
            handle.enable_write_to_file(f, u"radial_distribution", int(3))
            simulation.run(n_time_steps, 0.02)
            handle.flush()

        with h5py.File(fname, "r") as f2:
            bin_centers = f2["readdy/observables/radial_distribution/bin_centers"][:]
            distribution = f2["readdy/observables/radial_distribution/distribution"][:]
            for t in range(n_time_steps):
                np.testing.assert_equal(bin_centers, np.array(callback_centers[t]))
                np.testing.assert_equal(distribution[t], np.array(callback_rdf[t]))

    def test_histogram_along_axis_observable(self):
        fname = os.path.join(self.dir, "test_observables_hist_along_axis.h5")

        context = Context()
        context.kbt = 2
        context.pbc = [True, True, True]
        box_size = [10.,10.,10.]
        context.box_size = box_size
        context.particle_types.add("A", .2)
        context.particle_types.add("B", .2)
        context.potentials.add_harmonic_repulsion("A", "B", 10, 2.)
        simulation = Simulation("SingleCPU", context)

        simulation.add_particle("A", common.Vec(-2.5, 0, 0))
        simulation.add_particle("B", common.Vec(0, 0, 0))
        bin_borders = np.arange(0, 5, .01)
        n_time_steps = 50
        callback_hist = []

        def hist_callback(hist):
            callback_hist.append(hist)

        handle = simulation.register_observable_histogram_along_axis(2, bin_borders, 0, ["A", "B"], hist_callback)
        with closing(io.File.create(fname)) as f:
            handle.enable_write_to_file(f, u"hist_along_x_axis", int(3))
            simulation.run(n_time_steps, 0.02)
            handle.flush()

        with h5py.File(fname, "r") as f2:
            histogram = f2["readdy/observables/hist_along_x_axis/data"][:]
            time_series = f2["readdy/observables/hist_along_x_axis/time"]
            np.testing.assert_equal(time_series, np.array(range(0, n_time_steps+1))[::2])
            for t in range(n_time_steps // 2):
                np.testing.assert_equal(histogram[t], np.array(callback_hist[t]))

    def test_n_particles_observable(self):
        fname = os.path.join(self.dir, "test_observables_n_particles.h5")

        context = Context()
        box_size = [10.,10.,10.]
        context.kbt = 2
        context.pbc = [True, True, True]
        context.box_size = box_size
        context.particle_types.add("A", .2)
        context.particle_types.add("B", .2)
        simulation = Simulation("SingleCPU", context)

        simulation.add_particle("A", common.Vec(-2.5, 0, 0))
        simulation.add_particle("B", common.Vec(0, 0, 0))
        n_time_steps = 50
        callback_n_particles_a_b = []
        callback_n_particles_all = []

        def callback_ab(value):
            callback_n_particles_a_b.append(value)
            simulation.add_particle("A", common.Vec(-1, -1, -1))

        def callback_all(hist):
            callback_n_particles_all.append(hist)
            simulation.add_particle("A", common.Vec(-1, -1, -1))
            simulation.add_particle("B", common.Vec(-1, -1, -1))

        handle_a_b_particles = simulation.register_observable_n_particles(1, ["A", "B"], callback_ab)
        handle_all = simulation.register_observable_n_particles(1, [], callback_all)
        with closing(io.File.create(fname)) as f:
            handle_a_b_particles.enable_write_to_file(f, u"n_a_b_particles", int(3))
            handle_all.enable_write_to_file(f, u"n_particles", int(5))
            simulation.run(n_time_steps, 0.02)
            handle_all.flush()
            handle_a_b_particles.flush()

        with h5py.File(fname, "r") as f2:
            n_a_b_particles = f2["readdy/observables/n_a_b_particles/data"][:]
            n_particles = f2["readdy/observables/n_particles/data"][:]
            time_series = f2["readdy/observables/n_a_b_particles/time"]
            np.testing.assert_equal(time_series, np.array(range(0, n_time_steps+1)))
            for t in range(n_time_steps):
                np.testing.assert_equal(n_a_b_particles[t][0], callback_n_particles_a_b[t][0])
                np.testing.assert_equal(n_a_b_particles[t][1], callback_n_particles_a_b[t][1])
                np.testing.assert_equal(n_particles[t][0], callback_n_particles_all[t][0])

    def test_reactions_observable(self):
        fname = os.path.join(self.dir, "test_observables_particle_reactions.h5")
        context = Context()
        context.box_size = [10.,10.,10.]
        context.particle_types.add("A", .0)
        context.particle_types.add("B", .0)
        context.particle_types.add("C", .0)
        context.reactions.add_conversion("mylabel", "A", "B", .00001)
        context.reactions.add_conversion("A->B", "A", "B", 1.)
        context.reactions.add_fusion("B+C->A", "B", "C", "A", 1.0, 1.0, .5, .5)
        sim = Simulation("CPU", context)
        sim.add_particle("A", common.Vec(0, 0, 0))
        sim.add_particle("B", common.Vec(1.0, 1.0, 1.0))
        sim.add_particle("C", common.Vec(1.1, 1.0, 1.0))

        n_timesteps = 1

        handle = sim.register_observable_reactions(1)
        with closing(io.File.create(fname)) as f:
            handle.enable_write_to_file(f, u"reactions", int(3))
            loop = sim.create_loop(1)
            loop.write_config_to_file(f)
            loop.run(n_timesteps)

        type_str_to_id = {k: x["type_id"] for k, x in ioutils.get_particle_types(fname).items()}

        traj = readdy.Trajectory(fname)
        reac = traj.read_observable_reactions()
        with h5py.File(fname, "r") as f2:
            data = f2["readdy/observables/reactions"]
            time_series = f2["readdy/observables/reactions/time"]
            np.testing.assert_equal(time_series, np.array(range(0, n_timesteps+1)))

            def get_item(name, collection):
                return next(x for x in collection if x["name"] == name)

            import readdy.util.io_utils as io_utils
            reactions = io_utils.get_reactions(fname)

            mylabel_reaction = get_item(b"mylabel", reactions.values())
            np.testing.assert_allclose(mylabel_reaction["rate"], .00001)
            np.testing.assert_equal(mylabel_reaction["n_educts"], 1)
            np.testing.assert_equal(mylabel_reaction["n_products"], 1)
            np.testing.assert_equal(mylabel_reaction["educt_types"], [type_str_to_id["A"], 0])
            np.testing.assert_equal(mylabel_reaction["product_types"], [type_str_to_id["B"], 0])
            atob_reaction = get_item(b"A->B", reactions.values())
            np.testing.assert_equal(atob_reaction["rate"], 1.)
            np.testing.assert_equal(atob_reaction["n_educts"], 1)
            np.testing.assert_equal(atob_reaction["n_products"], 1)
            np.testing.assert_equal(mylabel_reaction["educt_types"], [type_str_to_id["A"], 0])
            np.testing.assert_equal(mylabel_reaction["product_types"], [type_str_to_id["B"], 0])

            fusion_reaction = get_item(b"B+C->A", reactions.values())
            np.testing.assert_equal(fusion_reaction["rate"], 1.)
            np.testing.assert_equal(fusion_reaction["educt_distance"], 1.)
            np.testing.assert_equal(fusion_reaction["n_educts"], 2)
            np.testing.assert_equal(fusion_reaction["n_products"], 1)
            np.testing.assert_equal(fusion_reaction["educt_types"], [type_str_to_id["B"], type_str_to_id["C"]])
            np.testing.assert_equal(fusion_reaction["product_types"], [type_str_to_id["A"], 0])

        _, records = reac
        np.testing.assert_equal(len(records), 2)
        # records of 1st time step
        for record in records[1]:
            np.testing.assert_(record.type in ('conversion', 'fusion'))
            if record.type == 'conversion':
                np.testing.assert_equal(record.position, np.array([.0, .0, .0]))
                np.testing.assert_equal(record.reaction_label, 'A->B')
            else:
                # fusion
                np.testing.assert_allclose(record.position, np.array([1.05, 1.0, 1.0]))
                np.testing.assert_equal(record.reaction_label, 'B+C->A')

    def test_reaction_counts_observable(self):
        fname = os.path.join(self.dir, "test_observables_particle_reaction_counts.h5")
        context = Context()
        context.box_size = [10., 10., 10.]
        context.particle_types.add("A", .0)
        context.particle_types.add("B", .0)
        context.particle_types.add("C", .0)
        context.reactions.add_conversion("mylabel", "A", "B", .00001)
        context.reactions.add_conversion("A->B", "A", "B", 1e16)
        context.reactions.add_fusion("B+C->A", "B", "C", "A", 1e16, 1.0, .5, .5)
        sim = Simulation("CPU", context)
        sim.add_particle("A", common.Vec(0, 0, 0))
        sim.add_particle("B", common.Vec(1.0, 1.0, 1.0))
        sim.add_particle("C", common.Vec(1.1, 1.0, 1.0))

        n_timesteps = 1
        handle = sim.register_observable_reaction_counts(1)
        with closing(io.File.create(fname)) as f:
            handle.enable_write_to_file(f, u"reactions", int(3))
            loop = sim.create_loop(1)
            loop.use_reaction_scheduler("Gillespie")
            loop.write_config_to_file(f)
            loop.run(1)

        import readdy.util.io_utils as io_utils
        reactions = io_utils.get_reactions(fname)

        with h5py.File(fname, "r") as f2:
            data = f2["readdy/observables/reactions"]
            time_series = f2["readdy/observables/reactions/time"]
            np.testing.assert_equal(time_series, np.array(range(0, n_timesteps+1)))

            def get_item(name, collection):
                return next(x for x in collection if x["name"] == name)

            mylabel_id = get_item(b"mylabel", reactions.values())["id"]
            atob_id = get_item(b"A->B", reactions.values())["id"]
            fusion_id = get_item(b"B+C->A", reactions.values())["id"]

            # counts of first time step, time is first index
            np.testing.assert_equal(data["counts/"+str(mylabel_id)][0], np.array([0]))
            np.testing.assert_equal(data["counts/"+str(atob_id)][0], np.array([0]))
            np.testing.assert_equal(data["counts/"+str(fusion_id)][0], np.array([0]))
            # counts of second time step
            np.testing.assert_equal(data["counts/"+str(mylabel_id)][1], np.array([0]))
            np.testing.assert_equal(data["counts/"+str(atob_id)][1], np.array([1]))
            np.testing.assert_equal(data["counts/"+str(fusion_id)][1], np.array([1]))

    def test_forces_observable(self):
        fname = os.path.join(self.dir, "test_observables_particle_forces.h5")
        context = Context()
        context.box_size = [13.,13.,13.]
        context.particle_types.add("A", .1)
        sim = Simulation("CPU", context)
        sim.add_particle("A", common.Vec(0, 0, 0))
        # every time step, add one particle
        sim.register_observable_n_particles(1, ["A"], lambda n: sim.add_particle("A", common.Vec(1.5, 2.5, 3.5)))
        handle = sim.register_observable_forces(1, [])
        n_timesteps = 19
        with closing(io.File.create(fname)) as f:
            handle.enable_write_to_file(f, u"forces", int(3))
            sim.run(n_timesteps, 1.)
            handle.flush()

        with h5py.File(fname, "r") as f2:
            data = f2["readdy/observables/forces/data"][:]
            time_series = f2["readdy/observables/forces/time"]
            np.testing.assert_equal(len(data), n_timesteps + 1)
            np.testing.assert_equal(time_series, np.array(range(0, n_timesteps+1)))
            for t, forces in enumerate(data):
                # we begin with two particles
                np.testing.assert_equal(len(forces), t + 2)
                np.testing.assert_equal(forces[0]["x"], 0)
                np.testing.assert_equal(forces[0]["y"], 0)
                np.testing.assert_equal(forces[0]["z"], 0)
                for i in range(1, len(forces)):
                    np.testing.assert_equal(forces[i]["x"], 0)
                    np.testing.assert_equal(forces[i]["y"], 0)
                    np.testing.assert_equal(forces[i]["z"], 0)


if __name__ == '__main__':
    unittest.main()
