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
Created on 19.6.17

@author: clonker
"""

from __future__ import print_function

import unittest

import numpy as np
import readdy._internal.readdybinding.api.top as top
import readdy._internal.readdybinding.common as common
from readdy._internal.readdybinding.api import KernelProvider
from readdy._internal.readdybinding.api import ParticleTypeFlavor
from readdy._internal.readdybinding.api import Simulation
from readdy.util import platform_utils
from readdy.util.testing_utils import ReaDDyTestCase

from readdy._internal.readdybinding.api import AnglePotentialConfiguration
from readdy._internal.readdybinding.api import BondedPotentialConfiguration
from readdy._internal.readdybinding.api import TorsionPotentialConfiguration

import readdy


class TestTopologyReactions(ReaDDyTestCase):
    @classmethod
    def setUpClass(cls):
        cls.kernel_provider = KernelProvider.get()
        cls.kernel_provider.load_from_dir(platform_utils.get_readdy_plugin_dir())

    def test_chain_decay_scpu(self):
        self.chain_decay("SingleCPU")

    def test_chain_decay_cpu(self):
        self.chain_decay("CPU")

    def _get_split_reaction(self):
        def reaction_function(topology):
            recipe = top.Recipe(topology)
            if topology.get_n_particles() > 1:
                edge = np.random.randint(0, topology.get_n_particles() - 1)
                recipe.remove_edge(edge, edge + 1)
            return recipe

        def rate_function(topology):
            if topology.get_n_particles() > 1:
                return topology.get_n_particles() / 20.
            else:
                return .0

        fun1 = top.ReactionFunction(reaction_function)
        fun2 = top.RateFunction(rate_function)

        reaction = top.StructuralTopologyReaction(fun1, fun2)
        reaction.roll_back_if_invalid()
        reaction.create_child_topologies_after_reaction()
        return reaction

    def _get_decay_reaction(self):
        def reaction_function(topology):
            recipe = top.Recipe(topology)
            if topology.get_n_particles() == 1:
                recipe.change_particle_type(0, "B")
            return recipe

        def rate_function(topology):
            return 1.0 if topology.get_n_particles() == 1 else 0

        fun1, fun2 = top.ReactionFunction(reaction_function), top.RateFunction(rate_function)
        reaction = top.StructuralTopologyReaction(fun1, fun2)
        reaction.raise_if_invalid()
        reaction.create_child_topologies_after_reaction()
        return reaction

    def chain_decay(self, kernel):
        sim = Simulation(kernel)
        sim.context.box_size = [10., 10., 10.]
        sim.context.topologies.add_type("TA")
        np.testing.assert_equal(sim.kernel_supports_topologies(), True)

        sim.context.particle_types.add("B", 1.0, ParticleTypeFlavor.NORMAL)
        sim.context.particle_types.add("Topology A", 1.0, ParticleTypeFlavor.TOPOLOGY)
        sim.context.topologies.configure_bond_potential("Topology A", "Topology A", BondedPotentialConfiguration(10, 10, "harmonic"))

        n_elements = 50.
        particles = [sim.create_topology_particle("Topology A", common.Vec(-5. + i * 10. / n_elements, 0, 0))
                     for i in range(int(n_elements))]
        topology = sim.add_topology("TA", particles)

        for i in range(int(n_elements - 1)):
            topology.get_graph().add_edge(i, i + 1)

        sim.context.topologies.add_structural_reaction("TA", self._get_decay_reaction())
        sim.context.topologies.add_structural_reaction("TA", self._get_split_reaction())

        # h = sim.register_observable_n_particles(1, [], lambda x: print("n particles=%s" % x))

        np.testing.assert_equal(1, len(sim.current_topologies()))

        sim.run(500, 1.)

        np.testing.assert_equal(0, len(sim.current_topologies()))

    def test_edges_decay(self):
        def dissociation_reaction_function(topology):
            recipe = readdy.StructuralReactionRecipe(topology)
            edges = topology.get_graph().get_edges()
            edge = edges[1 + np.random.randint(0, len(edges)-2)]
            recipe.remove_edge(edge)
            recipe.change_particle_type(edge[0].get().particle_index, "Head")
            recipe.change_particle_type(edge[1].get().particle_index, "Head")
            return recipe

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


if __name__ == '__main__':
    unittest.main()
