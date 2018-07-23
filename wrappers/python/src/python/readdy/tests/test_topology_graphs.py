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
Created on 24.03.17

@author: clonker
"""

from __future__ import print_function

import unittest

import numpy as np
import readdy._internal.readdybinding.common as common
from readdy._internal.readdybinding.api import BondedPotentialConfiguration
from readdy._internal.readdybinding.api import ParticleTypeFlavor
from readdy._internal.readdybinding.api import Simulation

from readdy.util.testing_utils import ReaDDyTestCase


class TestTopologyGraphs(ReaDDyTestCase):
    def test_sanity(self):
        sim = Simulation("SingleCPU")
        sim.context.box_size = [10., 10., 10.]
        np.testing.assert_equal(sim.kernel_supports_topologies(), True)
        sim.context.topologies.add_type("TA")
        sim.context.particle_types.add("T", 1.0, flavor=ParticleTypeFlavor.TOPOLOGY)
        sim.context.topologies.configure_bond_potential("T", "T", BondedPotentialConfiguration(10., 11., "harmonic"))
        particles = [sim.create_topology_particle("T", common.Vec(x, 0, 0)) for x in range(4)]
        top = sim.add_topology("TA", particles)
        graph = top.get_graph()
        graph.add_edge(0, 1)
        graph.add_edge(1, 2)
        graph.add_edge(2, 3)
        np.testing.assert_equal(len(graph.get_vertices()), 4)
        for v in graph.get_vertices():
            if v.particle_index == 0:
                np.testing.assert_equal(top.position_of_vertex(v), common.Vec(0, 0, 0))
                np.testing.assert_equal(len(v.neighbors()), 1)
                np.testing.assert_equal(1 in [vv.get().particle_index for vv in v], True)
            if v.particle_index == 1:
                np.testing.assert_equal(top.position_of_vertex(v), common.Vec(1, 0, 0))
                np.testing.assert_equal(len(v.neighbors()), 2)
                np.testing.assert_equal(0 in [vv.get().particle_index for vv in v], True)
                np.testing.assert_equal(2 in [vv.get().particle_index for vv in v], True)
            if v.particle_index == 2:
                np.testing.assert_equal(top.position_of_vertex(v), common.Vec(2, 0, 0))
                np.testing.assert_equal(len(v.neighbors()), 2)
                np.testing.assert_equal(1 in [vv.get().particle_index for vv in v], True)
                np.testing.assert_equal(3 in [vv.get().particle_index for vv in v], True)
            if v.particle_index == 3:
                np.testing.assert_equal(top.position_of_vertex(v), common.Vec(3, 0, 0))
                np.testing.assert_equal(len(v.neighbors()), 1)
                np.testing.assert_equal(2 in [vv.get().particle_index for vv in v], True)
        top.configure()
        sim.run(0, 1)

    def test_unconnected_graph(self):
        sim = Simulation("SingleCPU")
        sim.context.topologies.add_type("TA")
        sim.context.box_size = [10., 10., 10.]
        np.testing.assert_equal(sim.kernel_supports_topologies(), True)
        sim.context.particle_types.add("T", 1.0, flavor=ParticleTypeFlavor.TOPOLOGY)
        sim.context.topologies.configure_bond_potential("T", "T", BondedPotentialConfiguration(10, 11, "harmonic"))
        particles = [sim.create_topology_particle("T", common.Vec(0, 0, 0)) for _ in range(4)]
        top = sim.add_topology("TA", particles)
        graph = top.get_graph()
        graph.add_edge(0, 1)
        graph.add_edge(1, 2)
        with (np.testing.assert_raises(ValueError)):
            top.configure()

    def test_unbonded_edge(self):
        sim = Simulation("SingleCPU")
        sim.context.box_size = [10., 10., 10.]
        sim.context.topologies.add_type("TA")
        np.testing.assert_equal(sim.kernel_supports_topologies(), True)
        sim.context.particle_types.add("T", 1.0, flavor=ParticleTypeFlavor.TOPOLOGY)
        sim.context.particle_types.add("D", 1.0, flavor=ParticleTypeFlavor.TOPOLOGY)
        sim.context.topologies.configure_bond_potential("T", "T", BondedPotentialConfiguration(10., 11., "harmonic"))
        particles = [sim.create_topology_particle("T", common.Vec(0, 0, 0)) for _ in range(3)]
        particles.append(sim.create_topology_particle("D", common.Vec(0, 0, 0)))
        top = sim.add_topology("TA", particles)
        graph = top.get_graph()
        graph.add_edge(0, 1)
        graph.add_edge(1, 2)
        graph.add_edge(2, 3)
        with (np.testing.assert_raises(ValueError)):
            top.configure()

if __name__ == '__main__':
    unittest.main()
