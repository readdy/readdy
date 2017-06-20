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


class TestTopologyReactions(unittest.TestCase):
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
                return .00001

        fun1 = top.ReactionFunction(reaction_function)
        fun2 = top.RateFunction(rate_function)

        reaction = top.TopologyReaction(fun1, fun2)
        reaction.roll_back_if_invalid()
        reaction.create_child_topologies_after_reaction()
        return reaction

    def _get_decay_reaction(self, typeidb):
        def reaction_function(topology):
            recipe = top.Recipe(topology)
            if topology.get_n_particles() == 1:
                recipe.change_particle_type(0, typeidb)
            return recipe

        def rate_function(topology):
            return 1.0 if topology.get_n_particles() == 1 else 0

        fun1, fun2 = top.ReactionFunction(reaction_function), top.RateFunction(rate_function)
        reaction = top.TopologyReaction(fun1, fun2)
        reaction.raise_if_invalid()
        reaction.create_child_topologies_after_reaction()
        return reaction

    def chain_decay(self, kernel):
        sim = Simulation()
        sim.set_kernel(kernel)
        sim.box_size = common.Vec(10, 10, 10)
        np.testing.assert_equal(sim.kernel_supports_topologies(), True)

        typeid_b = sim.register_particle_type("B", 1.0, 1.0, ParticleTypeFlavor.NORMAL)
        sim.register_particle_type("Topology A", 1.0, 1.0, ParticleTypeFlavor.TOPOLOGY)
        sim.configure_topology_bond_potential("Topology A", "Topology A", 10, 10)

        n_elements = 50.
        particles = [sim.create_topology_particle("Topology A", common.Vec(-5. + i * 10. / n_elements, 0, 0))
                     for i in range(int(n_elements))]
        topology = sim.add_topology(particles)

        for i in range(int(n_elements - 1)):
            topology.get_graph().add_edge(i, i + 1)

        topology.add_reaction(self._get_decay_reaction(typeid_b))
        topology.add_reaction(self._get_split_reaction())

        # h = sim.register_observable_n_particles(1, [], lambda x: print("n particles=%s" % x))

        np.testing.assert_equal(1, len(sim.current_topologies()))

        sim.run_scheme_readdy(True).evaluate_topology_reactions().configure_and_run(int(500), float(1.0))

        np.testing.assert_equal(0, len(sim.current_topologies()))


if __name__ == '__main__':
    unittest.main()
