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
Created on 19.6.17

@author: clonker
"""

from __future__ import print_function

import unittest

import numpy as np
import readdy
import readdy._internal.readdybinding.api.top as top
import readdy._internal.readdybinding.common as common
from readdy._internal.readdybinding.api import KernelProvider
from readdy._internal.readdybinding.api import ParticleTypeFlavor
from readdy._internal.readdybinding.api import Simulation
from readdy._internal.readdybinding.api import Context
from readdy.util import platform_utils
from readdy.util.testing_utils import ReaDDyTestCase

from readdy._internal.readdybinding.api import AnglePotentialConfiguration
from readdy._internal.readdybinding.api import BondedPotentialConfiguration
from readdy._internal.readdybinding.api import TorsionPotentialConfiguration


class TestTopologyReactions(ReaDDyTestCase):

    def test_chain_decay_scpu(self):
        self.chain_decay("SingleCPU")

    def test_chain_decay_cpu(self):
        self.chain_decay("CPU")

    def _get_split_reaction(self):
        def reaction_function(topology):
            recipe = readdy.StructuralReactionRecipe(topology)
            if topology.get_n_particles() > 1:
                edge = np.random.randint(0, topology.get_n_particles() - 1)
                recipe.remove_edge(edge, edge + 1)
            return recipe._get()

        def rate_function(topology):
            if topology.get_n_particles() > 1:
                return topology.get_n_particles() / 20.
            else:
                return .0

        fun1 = top.ReactionFunction(reaction_function)
        fun2 = top.RateFunction(rate_function)

        reaction = top.StructuralTopologyReaction("split", fun1, fun2)
        reaction.create_child_topologies_after_reaction()
        return reaction

    def _get_decay_reaction(self):
        def reaction_function(topology):
            recipe = readdy.StructuralReactionRecipe(topology)
            if topology.get_n_particles() == 1:
                recipe.change_particle_type(0, "B")
            return recipe._get()

        def rate_function(topology):
            return 1.0 if topology.get_n_particles() == 1 else 0

        fun1, fun2 = top.ReactionFunction(reaction_function), top.RateFunction(rate_function)
        reaction = top.StructuralTopologyReaction("decay", fun1, fun2)
        reaction.create_child_topologies_after_reaction()
        return reaction

    def chain_decay(self, kernel):
        context = Context()
        context.box_size = [10., 10., 10.]
        context.topologies.add_type("TA")

        context.particle_types.add("B", 1.0, ParticleTypeFlavor.NORMAL)
        context.particle_types.add("Topology A", 1.0, ParticleTypeFlavor.TOPOLOGY)
        context.topologies.configure_bond_potential("Topology A", "Topology A", BondedPotentialConfiguration(10, 10, "harmonic"))

        context.topologies.add_structural_reaction("TA", self._get_decay_reaction())
        context.topologies.add_structural_reaction("TA", self._get_split_reaction())

        sim = Simulation(kernel, context)
        np.testing.assert_equal(sim.kernel_supports_topologies(), True)

        n_elements = 50.
        particles = [sim.create_topology_particle("Topology A", common.Vec(-5. + i * 10. / n_elements, 0, 0))
                     for i in range(int(n_elements))]
        topology = sim.add_topology("TA", particles)

        for i in range(int(n_elements - 1)):
            topology.get_graph().add_edge(i, i + 1)



        # h = sim.register_observable_n_particles(1, [], lambda x: print("n particles=%s" % x))

        np.testing.assert_equal(1, len(sim.current_topologies))

        sim.run(500, 1.)

        np.testing.assert_equal(0, len(sim.current_topologies))

    def test_edges_decay(self):
        def dissociation_reaction_function(topology):
            recipe = readdy.StructuralReactionRecipe(topology)
            edges = topology.get_graph().get_edges()
            edge = edges[1 + np.random.randint(0, len(edges)-2)]
            recipe.remove_edge(edge[0], edge[1])
            recipe.change_particle_type(edge[0].get().particle_index, "Head")
            recipe.change_particle_type(edge[1].get().particle_index, "Head")
            return recipe._get()

        system = readdy.ReactionDiffusionSystem(box_size=[100, 100, 100])
        system.topologies.add_type("Polymer")
        system.add_topology_species("Head", .002)
        system.add_topology_species("Tail", .002)
        simulation = system.simulation(kernel="SingleCPU")
        head_position = [0, 0, 0]
        tail1 = [0, 0, 1]
        tail2 = [0, 0, 2]
        head_position2 = [0, 0, 3]
        top = simulation.add_topology("Polymer", ["Head", "Tail", "Tail", "Head"],
                                      np.array([head_position, tail1, tail2, head_position2]).squeeze())
        top.get_graph().add_edge(0, 1)
        top.get_graph().add_edge(1, 2)
        top.get_graph().add_edge(2, 3)

        dissociation_reaction_function(top)

    def test_spatial_reaction_rate_function_SingleCPU(self):
        self.spatial_reaction_rate_function("SingleCPU")

    def test_spatial_reaction_rate_function_CPU(self):
        self.spatial_reaction_rate_function("CPU")

    @staticmethod
    def spatial_reaction_rate_function(kernel):
        """
        Create a small simulation of 4 polymers that
        form arranged as a square. Heads of the polymers
        are close, such that they could bind.
        The rate function should allow only two
        of the polymers to bind, so that we end up with
        2 pairs of two connected polymers.
        """
        def rate_function(top1, top2):
            vert1 = top1.get_graph().get_vertices()
            vert2 = top2.get_graph().get_vertices()
            if len(vert1) + len(vert2) > 12:
                return 0.0
            return 1e10

        system = readdy.ReactionDiffusionSystem(box_size=[30., 30., 30.])
        system.topologies.add_type("Polymer")
        system.add_topology_species("Head", 0.002)
        system.add_topology_species("Core", 0.002)

        system.topologies.configure_harmonic_bond("Head", "Core", force_constant=50, length=1.)
        system.topologies.configure_harmonic_bond("Core", "Core", force_constant=50, length=1.)

        system.topologies.add_spatial_reaction(
            "Association: Polymer(Head) + Polymer(Head) -> Polymer(Core--Core)",
            rate=rate_function, radius=2.0
        )

        simulation = system.simulation(kernel=kernel)
        types_and_positions = TestTopologyReactions._get_polymer_types_and_positions()
        for t, p in types_and_positions:
            top = simulation.add_topology("Polymer", t, p)
            for i in range(5):
                top.get_graph().add_edge(i, i+1)

        simulation.run(10, 1.)

        np.testing.assert_equal(2, len(simulation.current_topologies))

    @staticmethod
    def _get_polymer_types_and_positions():
        """
        Construct a square of 4 polymers.
        """
        types_and_positions = [TestTopologyReactions._get_types_and_positions_polymer_1(),
                               TestTopologyReactions._get_types_and_positions_polymer_2(),
                               TestTopologyReactions._get_types_and_positions_polymer_3(),
                               TestTopologyReactions._get_types_and_positions_polymer_4()]

        return types_and_positions

    @staticmethod
    def _get_types_and_positions_polymer_1():
        types = ["Head"] + ["Core"]*4 + ["Head"]
        positions = np.zeros((6, 3))
        positions[:, 0] = np.arange(6)
        return types, positions

    @staticmethod
    def _get_types_and_positions_polymer_2():
        types = ["Head"] + ["Core"]*4 + ["Head"]
        positions = np.zeros((6, 3))
        positions[:, 0] = 5
        positions[:, 1] = np.arange(0, -6, -1)
        positions[:, 2] = 1
        return types,positions

    @staticmethod
    def _get_types_and_positions_polymer_3():
        types = ["Head"] + ["Core"]*4 + ["Head"]
        positions = np.zeros((6, 3))
        positions[:, 0] = np.arange(5, -1, -1)
        positions[:, 1] = -5
        return types, positions

    @staticmethod
    def _get_types_and_positions_polymer_4():
        types = ["Head"] + ["Core"]*4 + ["Head"]
        positions = np.zeros((6, 3))
        positions[:, 1] = np.arange(6) - 5
        positions[:, 2] = 1
        return types, positions



if __name__ == '__main__':
    unittest.main()
