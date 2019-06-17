import os
import shutil

import readdy
import tempfile
import unittest

import numpy as np

class TestTopologyReactionCount(unittest.TestCase):

    @classmethod
    def setUpClass(cls) -> None:
        cls.dir = tempfile.mkdtemp("test-topology-reaction-count")

    @classmethod
    def tearDownClass(cls) -> None:
        shutil.rmtree(cls.dir, ignore_errors=True)

    def _test_kernel(self, kernel):

        system = readdy.ReactionDiffusionSystem(box_size=[10, 10, 10])
        system.topologies.add_type("T1")
        system.topologies.add_type("T2")
        system.add_species("A")
        system.add_topology_species("B")
        system.topologies.configure_harmonic_bond("B", "B", 1., .1)
        system.add_topology_species("C")
        system.topologies.configure_harmonic_bond("C", "C", 1., .1)
        system.topologies.add_spatial_reaction("attach: T1(B) + (A) -> T1(B--B)", rate=1e-1, radius=.5)

        def flip1(topology):
            recipe = readdy.StructuralReactionRecipe(topology)
            for v in topology.graph.vertices:
                recipe.change_particle_type(v, "C")
            recipe.change_topology_type("T2")
            return recipe

        def flip2(topology):
            recipe = readdy.StructuralReactionRecipe(topology)
            for v in topology.graph.vertices:
                recipe.change_particle_type(v, "B")
            recipe.change_topology_type("T1")
            return recipe

        system.topologies.add_structural_reaction("flip_types_1", "T1", flip1, lambda x: 5e-2)
        system.topologies.add_structural_reaction("flip_types_2", "T2", flip2, lambda x: 5e-2)

        sim = system.simulation(kernel=kernel)
        sim.output_file = os.path.join(self.dir, "out_{}.h5".format(kernel))

        collected_counts = []

        def callback(results):
            nonlocal collected_counts
            collected_counts.append(results)

        sim.observe.reaction_counts(1, callback=callback)
        sim.observe.number_of_particles(1, types=["A", "B", "C"])

        sim.add_particles("A", np.random.normal(scale=1, size=(1000, 3)))

        for _ in range(10):
            sim.add_topology("T1", "B", np.random.normal(size=(1, 3)))

        sim.run(1000, timestep=1, show_summary=False)

        traj = readdy.Trajectory(sim.output_file)

        times, n_particles = traj.read_observable_number_of_particles()
        times2, counts = traj.read_observable_reaction_counts()

        np.testing.assert_array_equal(times, times2)

        assert not counts["reactions"]

        for t, (cA, cB, cC), cc in zip(times, n_particles, collected_counts):
            np.testing.assert_equal(cA + cB + cC, 1010)
            cc_normal = cc[0]
            assert not cc_normal
            cc_spatial = cc[1]
            cc_structural = cc[2]

            for sp in cc_spatial.keys():
                recorded = counts["spatial_topology_reactions"][sp][t]
                assert cc_spatial[sp] == recorded, f"Got {cc_spatial[sp]} != {recorded} (t={t})"
            for st in cc_structural.keys():
                recorded = counts["structural_topology_reactions"][st][t]
                assert cc_structural[st] == recorded, f"Got {cc_structural[st]} != {recorded} (t={t})"

    def test_scpu(self):
        self._test_kernel("SingleCPU")

    def test_cpu(self):
        self._test_kernel("CPU")


if __name__ == '__main__':
    unittest.main()
