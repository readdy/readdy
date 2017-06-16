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
import readdy._internal.readdybinding.common as common
import readdy._internal.readdybinding.api.top as top
from readdy._internal.readdybinding.api import KernelProvider
from readdy._internal.readdybinding.api import ParticleTypeFlavor
from readdy._internal.readdybinding.api import Simulation
from readdy.util import platform_utils

class ReactionGenerator(top.ReactionFunctionGenerator):

    def __init__(self):
        super(ReactionGenerator, self).__init__()

    def generate(self):
        print("foooooooo")

        def reaction_function():
            print("foo")
            recipe = top.Recipe(topology)
            if topology.get_n_particles() > 1:
                edge = np.random.randint(0, topology.get_n_particles()-2)
                recipe.remove_edge(edge, edge+1)
            return recipe
        return reaction_function

class RateGenerator(top.RateFunctionGenerator):

    def __init__(self):
        super(RateGenerator, self).__init__()

    def generate(self, topology):
        print("aga!")

        def rate_function():
            print("huhu")
            if topology.get_n_particles() > 1:
                return topology.get_n_particles() / 50.
            else:
                return 0
        return rate_function


class TestTopologyReactions(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.kernel_provider = KernelProvider.get()
        cls.kernel_provider.load_from_dir(platform_utils.get_readdy_plugin_dir())

    def test_chain_decay_scpu(self):
        self.chain_decay("SingleCPU")

    def test_chain_decay_cpu(self):
        pass # self.chain_decay("CPU")

    def _get_decay_reaction(self):





        reaction = top.TopologyReaction.create(ReactionGenerator(), RateGenerator())
        reaction.roll_back_if_invalid()
        reaction.create_child_topologies_after_reaction()
        return reaction

    def chain_decay(self, kernel):
        sim = Simulation()
        sim.set_kernel(kernel)
        sim.box_size = common.Vec(10, 10, 10)
        np.testing.assert_equal(sim.kernel_supports_topologies(), True)

        sim.register_particle_type("Topology A", 1.0, 1.0, ParticleTypeFlavor.TOPOLOGY)
        sim.register_particle_type("B", 1.0, 1.0, ParticleTypeFlavor.NORMAL)
        sim.configure_topology_bond_potential("Topology A", "Topology A", 10, 10)

        n_elements = 50.
        particles = [sim.create_topology_particle("Topology A", common.Vec(-5. + i * 10. / n_elements, 0, 0))
                     for i in range(int(n_elements))]
        topology = sim.add_topology(particles)

        for i in range(int(n_elements-1)):
            topology.get_graph().add_edge(i, i+1)

        topology.add_reaction(self._get_decay_reaction())

        sim.run_scheme_readdy(True).evaluate_topology_reactions().configure_and_run(int(5), float(1.0))


if __name__ == '__main__':
    unittest.main()
