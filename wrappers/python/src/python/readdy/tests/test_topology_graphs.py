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
Created on 24.03.17

@author: clonker
"""

from __future__ import print_function

import unittest
import numpy as np
import readdy._internal.readdybinding.common as common
from readdy._internal.readdybinding.api import Simulation
from readdy._internal.readdybinding.api import ParticleTypeFlavor


class TestTopologyGraphs(unittest.TestCase):
    def test_sanity(self):
        sim = Simulation()
        sim.set_kernel("SingleCPU")
        sim.box_size = common.Vec(10, 10, 10)
        np.testing.assert_equal(sim.kernel_supports_topologies(), True)
        sim.register_particle_type("T", 1.0, .5, flavor=ParticleTypeFlavor.TOPOLOGY)
        sim.configure_topology_bond_potential("T", "T", 10., 11.)
        particles = [sim.create_topology_particle("T", common.Vec(0, 0, 0)) for _ in range(4)]
        labels = ["%s" % i for i in range(4)]
        top = sim.add_topology(particles, labels)
        graph = top.get_graph()
        graph.add_edge("0", "1")
        graph.add_edge(1, 2)
        graph.add_edge("2", "3")
        np.testing.assert_equal(len(graph.get_vertices()), 4)
        for v in graph.get_vertices():
            if v.label == "0":
                np.testing.assert_equal(len(v.neighbors()), 1)
                np.testing.assert_equal(1 in [vv.get().particle_index for vv in v], True)
            if v.label == "1":
                np.testing.assert_equal(len(v.neighbors()), 2)
                np.testing.assert_equal(0 in [vv.get().particle_index for vv in v], True)
                np.testing.assert_equal(2 in [vv.get().particle_index for vv in v], True)
            if v.label == "2":
                np.testing.assert_equal(len(v.neighbors()), 2)
                np.testing.assert_equal(1 in [vv.get().particle_index for vv in v], True)
                np.testing.assert_equal(3 in [vv.get().particle_index for vv in v], True)
            if v.label == "3":
                np.testing.assert_equal(len(v.neighbors()), 1)
                np.testing.assert_equal(2 in [vv.get().particle_index for vv in v], True)
        top.configure()
        sim.run_scheme_readdy(True).configure_and_run(1, 0)

    def test_unconnected_graph(self):
        sim = Simulation()
        sim.set_kernel("SingleCPU")
        sim.box_size = common.Vec(10, 10, 10)
        np.testing.assert_equal(sim.kernel_supports_topologies(), True)
        sim.register_particle_type("T", 1.0, .5, flavor=ParticleTypeFlavor.TOPOLOGY)
        sim.configure_topology_bond_potential("T", "T", 10., 11.)
        particles = [sim.create_topology_particle("T", common.Vec(0, 0, 0)) for _ in range(4)]
        top = sim.add_topology(particles)
        graph = top.get_graph()
        graph.add_edge(0, 1)
        graph.add_edge(1, 2)
        with (np.testing.assert_raises(ValueError)):
            top.configure()

    def test_unbonded_edge(self):
        sim = Simulation()
        sim.set_kernel("SingleCPU")
        sim.box_size = common.Vec(10, 10, 10)
        np.testing.assert_equal(sim.kernel_supports_topologies(), True)
        sim.register_particle_type("T", 1.0, .5, flavor=ParticleTypeFlavor.TOPOLOGY)
        sim.register_particle_type("D", 1.0, .5, flavor=ParticleTypeFlavor.TOPOLOGY)
        sim.configure_topology_bond_potential("T", "T", 10., 11.)
        particles = [sim.create_topology_particle("T", common.Vec(0, 0, 0)) for _ in range(3)]
        particles.append(sim.create_topology_particle("D", common.Vec(0, 0, 0)))
        top = sim.add_topology(particles)
        graph = top.get_graph()
        graph.add_edge(0, 1)
        graph.add_edge(1, 2)
        graph.add_edge(2, 3)
        with (np.testing.assert_raises(ValueError)):
            top.configure()

if __name__ == '__main__':
    unittest.main()
