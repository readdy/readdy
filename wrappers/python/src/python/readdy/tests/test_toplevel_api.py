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
Created on 26.09.17

@author: clonker
"""
import shutil
import tempfile
import os

import numpy as np

import readdy
from readdy.util.testing_utils import ReaDDyTestCase


class TestTopLevelAPI(ReaDDyTestCase):
    """
    A bunch of sanity checks
    """

    def test_kbt(self):
        rdf = readdy.ReactionDiffusionSystem()
        rdf.kbt = 5.
        np.testing.assert_equal(rdf.kbt, 5.)

    def test_box_size(self):
        rdf = readdy.ReactionDiffusionSystem()
        rdf.box_size = [1., 2., 3.]
        np.testing.assert_equal(rdf.box_size, [1., 2., 3.])
        rdf.box_size = np.array([5., 6., 7.])
        np.testing.assert_equal(rdf.box_size, [5., 6., 7.])
        rdf.box_size = (1., 5., 7.)
        np.testing.assert_equal(rdf.box_size, [1., 5., 7.])
        np.testing.assert_equal(rdf.box_volume, 5. * 7.)

    def test_pbc(self):
        rdf = readdy.ReactionDiffusionSystem()
        rdf.periodic_boundary_conditions = True, False, True
        np.testing.assert_equal(rdf.periodic_boundary_conditions, [True, False, True])
        rdf.periodic_boundary_conditions = np.array([False, False, True])
        np.testing.assert_equal(rdf.periodic_boundary_conditions, [False, False, True])

    def test_species(self):
        rdf = readdy.ReactionDiffusionSystem()
        rdf.add_species("A", 1.)
        self.assertTrue("A" in rdf.registered_species())
        rdf.add_topology_species("Top A", 10.)
        self.assertTrue("A" in rdf.registered_species() and "Top A" in rdf.registered_species())

    def test_topology_potentials(self):
        rdf = readdy.ReactionDiffusionSystem()
        rdf.add_topology_species("A", 1.)
        rdf.add_topology_species("B", 1.)
        rdf.add_topology_species("C", 1.)
        rdf.add_topology_species("D", 1.)
        rdf.topologies.configure_harmonic_bond("A", "B", 1., 0.)
        rdf.topologies.configure_harmonic_angle("A", "B", "C", 1., 0.)
        rdf.topologies.configure_cosine_dihedral("A", "B", "C", "D", 1., 1, 0.)

    def test_spatial_topology_reactions(self):
        rdf = readdy.ReactionDiffusionSystem()
        rdf.add_topology_species("A", 0.)
        rdf.topologies.add_type("T1")
        rdf.topologies.add_type("T2")
        rdf.topologies.add_type("T3")
        rdf.topologies.add_spatial_reaction("test_fusion: T1(A)+T2(A) -> T3(A--A)", 1., 1.)
        rdf.topologies.add_spatial_reaction("test_enzymatic: T1(A)+T2(A) -> T3(A)+T2(A)", 1., 1.)

    def test_structural_topology_reactions(self):
        rdf = readdy.ReactionDiffusionSystem()
        rdf.add_topology_species("foo")
        rdf.topologies.add_type("foofoo")

        def reaction_fun(topology):
            return readdy.StructuralReactionRecipe(topology) \
                .change_particle_type(0, "foo")

        def rate_fun(topology):
            return len(topology.particles)

        rdf.topologies.add_structural_reaction("foofoo", reaction_fun, rate_fun)

    def test_potentials(self):
        rdf = readdy.ReactionDiffusionSystem()
        rdf.add_species("A")
        rdf.potentials.add_box("A", 1.0, [1.0, 1.0, 1.0], [1.0, 1.0, 1.0])
        rdf.potentials.add_harmonic_repulsion("A", "A", 1.0, 1.0)
        rdf.potentials.add_lennard_jones("A", "A", 12, 6, 10, True, 1, 1)
        rdf.potentials.add_screened_electrostatics("A", "A", 10, 10, 10, 10, 10, 10)
        rdf.potentials.add_sphere_in("A", 10, (10, 10, 10), 1)
        rdf.potentials.add_sphere_out("A", 10, (10, 10, 10), 1)
        rdf.potentials.add_spherical_barrier("A", 1, 1, (0, 0, 0), 1)
        rdf.potentials.add_weak_interaction_piecewise_harmonic("A", "A", 10, 10, 10, 10)

    def test_simulation(self):
        rdf = readdy.ReactionDiffusionSystem()
        rdf.add_species("A")
        _ = rdf.simulation("CPU")

    def test_observables(self):
        rdf = readdy.ReactionDiffusionSystem()
        rdf.add_species("A")
        simulation = rdf.simulation("CPU")
        simulation.observe.rdf(5, [0., 1., 2.], ["A"], "A", 10)
        simulation.observe.forces(5)
        simulation.observe.number_of_particles(5)
        simulation.observe.particle_positions(5)
        simulation.observe.particles(5)
        simulation.observe.reaction_counts(5)
        simulation.observe.reactions(5)
        simulation.run(10, .1)

    def test_add_particles(self):
        rdf = readdy.ReactionDiffusionSystem()
        rdf.box_size = [10., 10., 10.]
        rdf.add_species("A")
        sim = rdf.simulation("CPU")
        sim.add_particles("A", np.random.random((10000, 3)))

    def test_add_topology(self):
        rdf = readdy.ReactionDiffusionSystem()
        rdf.box_size = [10., 10., 10.]
        rdf.topologies.add_type("toptype")
        rdf.add_topology_species("TopA")
        rdf.add_topology_species("TopB")
        sim = rdf.simulation(kernel="SingleCPU")
        top1positions = np.random.random((4, 3))
        topology1 = sim.add_topology("toptype", "TopA", top1positions)
        for i, v in enumerate(topology1.get_graph().get_vertices()):
            np.testing.assert_equal("TopA", topology1.particle_type_of_vertex(v))
            np.testing.assert_equal(readdy.api.utils.vec3_of(top1positions[i, :]), topology1.position_of_vertex(v))
        top2_types = ["TopB"] + ["TopA" for _ in range(9)]
        top2positions = np.random.random((10, 3))
        topology2 = sim.add_topology("toptype", top2_types, top2positions)
        for i, v in enumerate(topology2.get_graph().get_vertices()):
            np.testing.assert_equal(readdy.api.utils.vec3_of(top2positions[i, :]), topology2.position_of_vertex(v))
            if i == 0:
                np.testing.assert_equal("TopB", topology2.particle_type_of_vertex(v))
            else:
                np.testing.assert_equal("TopA", topology2.particle_type_of_vertex(v))


class TestTopLevelAPIObservables(ReaDDyTestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.tempdir = tempfile.mkdtemp("top-level-test-obs")

    @classmethod
    def tearDownClass(cls):
        super().tearDownClass()
        shutil.rmtree(cls.tempdir, ignore_errors=True)

    def test_write_traj(self):
        traj_fname = os.path.join(self.tempdir, "traj.h5")

        rdf = readdy.ReactionDiffusionSystem()
        rdf.box_size = (10, 10, 10)
        rdf.add_species("A", diffusion_constant=1.0)
        rdf.reactions.add_conversion("myconversion", "A", "A", 1.0)
        rdf.reactions.add_fusion("myfusion", "A", "A", "A", 2, .5)
        rdf.potentials.add_harmonic_repulsion("A", "A", 1., .2)
        sim = rdf.simulation(kernel="SingleCPU")
        sim.output_file = traj_fname
        sim.record_trajectory(1)
        sim.add_particles("A", np.random.random((100, 3)))
        recorded_positions = []
        sim.observe.particle_positions(1, callback=lambda x: recorded_positions.append(x))
        sim.run(50, 1e-3)

        traj = readdy.Trajectory(traj_fname)

        np.testing.assert_equal(traj.diffusion_constants["A"], 1.0)
        np.testing.assert_("A" in traj.particle_types.keys())
        np.testing.assert_equal(len(traj.reactions), 2)

        conv = next(x for x in traj.reactions if x.name == "myconversion")
        np.testing.assert_equal(conv.type, "conversion")
        np.testing.assert_equal(conv.name, "myconversion")
        np.testing.assert_equal(conv.rate, 1.0)
        np.testing.assert_equal(conv.educt_types, ["A"])
        np.testing.assert_equal(conv.product_types, ["A"])

        fusion = next(x for x in traj.reactions if x.name == "myfusion")
        np.testing.assert_equal(fusion.type, "fusion")
        np.testing.assert_equal(fusion.name, "myfusion")
        np.testing.assert_equal(fusion.rate, 2)
        np.testing.assert_equal(fusion.educt_distance, .5)
        np.testing.assert_equal(fusion.educt_types, ["A", "A"])
        np.testing.assert_equal(fusion.product_types, ["A"])

        for idx, frame in enumerate(traj.read()):
            recorded = recorded_positions[idx]
            np.testing.assert_equal(len(recorded), len(frame))
            for e_idx, entry in enumerate(frame):
                pos = recorded[e_idx]
                np.testing.assert_equal(pos.toarray(), entry.position)
                np.testing.assert_equal("NORMAL", entry.flavor)
                np.testing.assert_equal("A", entry.type)
                np.testing.assert_equal(idx, entry.t)

    def run_readwrite_test_for(self, kernel, reaction_handler):
        traj_fname = os.path.join(self.tempdir, "traj_{}_{}.h5".format(kernel, reaction_handler))
        traj_fname2 = os.path.join(self.tempdir, "traj2_{}_{}.h5".format(kernel, reaction_handler))

        rdf = readdy.ReactionDiffusionSystem()
        rdf.box_size = (10, 10, 10)
        rdf.add_species("A", diffusion_constant=1.0)
        rdf.add_species("B", diffusion_constant=1.0)
        rdf.reactions.add_conversion("myconversion", "A", "B", 1.0)
        rdf.reactions.add_fusion("myfusion", "A", "A", "A", 2, .5)
        rdf.reactions.add_fission("myfission", "A", "A", "A", 2, .5)
        rdf.potentials.add_harmonic_repulsion("A", "A", 1., .2)
        sim = rdf.simulation(kernel=kernel, reaction_handler=reaction_handler)
        sim.output_file = traj_fname
        sim.add_particles("A", np.random.random((100, 3)))

        sim.observe.particle_positions(1)
        sim.observe.particles(1)
        sim.observe.rdf(1, bin_borders=np.arange(-5, 5, 1.), types_count_from=["A"], types_count_to=["A"],
                        particle_to_density=1. / rdf.box_volume)
        sim.observe.number_of_particles(1, types=["B", "A"])
        reactions = []
        sim.observe.reactions(1, callback=lambda x: reactions.append(x))
        sim.observe.reaction_counts(1)
        sim.observe.forces(1)
        sim.run(50, 1e-3)

        sim.output_file = traj_fname2
        sim.run(50, 1e-3)

        for fname in [traj_fname, traj_fname2]:
            traj = readdy.Trajectory(fname)

            time, positions = traj.read_observable_particle_positions()
            np.testing.assert_equal(len(time), 51)
            np.testing.assert_equal(len(positions), 51)

            time, types, ids, positions = traj.read_observable_particles()
            np.testing.assert_equal(len(time), 51)
            np.testing.assert_equal(len(types), 51)
            np.testing.assert_equal(len(ids), 51)
            np.testing.assert_equal(len(positions), 51)

            time, bin_centers, rdf = traj.read_observable_rdf()
            np.testing.assert_equal(len(time), 51)
            np.testing.assert_equal(len(bin_centers), len(np.arange(-5, 5, 1.)) - 1)
            np.testing.assert_equal(rdf.shape, (51, len(np.arange(-5, 5, 1)) - 1))

            time, counts = traj.read_observable_number_of_particles()
            np.testing.assert_equal(len(time), 51)
            np.testing.assert_equal(len(counts), 51)

            time, records = traj.read_observable_reactions()
            np.testing.assert_equal(len(time), 51)
            np.testing.assert_equal(len(records), 51)

            time, forces = traj.read_observable_forces()
            np.testing.assert_equal(len(time), 51)
            np.testing.assert_equal(len(forces), 51)

            time, counts = traj.read_observable_reaction_counts()
            np.testing.assert_equal(len(time), 51)
            np.testing.assert_equal(len(counts.keys()), 3)
            for t, rr, counts_1, counts_2, counts_3 in zip(time, records, counts["myconversion"],
                                                           counts["myfusion"], counts["myfission"]):
                convrecords = [r for r in rr if r.reaction_label == "myconversion"]
                fusrecords = [r for r in rr if r.reaction_label == "myfusion"]
                fissrecords = [r for r in rr if r.reaction_label == "myfission"]
                np.testing.assert_equal(counts_1, len(convrecords),
                                        err_msg="conversion count mismatch: t={}, counts={}, nrecords={}, kernel={}"
                                        .format(t, counts_1, len(convrecords), kernel))
                np.testing.assert_equal(counts_2, len(fusrecords),
                                        err_msg="fusion count mismatch: t={}, counts={}, nrecords={}, kernel={}"
                                        .format(t, counts_2, len(fusrecords), kernel))
                np.testing.assert_equal(counts_3, len(fissrecords),
                                        err_msg="fission count mismatch: t={}, counts={}, nrecords={}, kernel={}"
                                        .format(t, counts_3, len(fissrecords), kernel))

            for curr_positions, curr_types, curr_ids, curr_forces in zip(positions, types, ids, forces):
                np.testing.assert_equal(len(curr_positions), len(curr_types))
                np.testing.assert_equal(len(curr_types), len(curr_ids))
                np.testing.assert_equal(len(curr_ids), len(curr_forces))

    def test_readwrite_observables_singlecpu_uncontrolled(self):
        self.run_readwrite_test_for("SingleCPU", "UncontrolledApproximation")

    def test_readwrite_observables_singlecpu_gillespie(self):
        self.run_readwrite_test_for("SingleCPU", "Gillespie")

    def test_readwrite_observables_cpu_uncontrolled(self):
        self.run_readwrite_test_for("CPU", "UncontrolledApproximation")

    def test_readwrite_observables_cpu_gillespie(self):
        self.run_readwrite_test_for("CPU", "Gillespie")
