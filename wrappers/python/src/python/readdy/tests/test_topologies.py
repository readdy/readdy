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

from __future__ import print_function

"""
Created on 09.02.17

@author: clonker
"""
import unittest
import numpy as np
import readdy._internal.readdybinding.common as common
from readdy._internal.readdybinding.api import Simulation
from readdy._internal.readdybinding.api import ParticleTypeFlavor
from readdy._internal.readdybinding.api.top import HarmonicBondPotentialBond


class TestTopologies(unittest.TestCase):
    def test_sanity(self):
        sim = Simulation()
        sim.set_kernel("SingleCPU")
        sim.register_topology_type("TA")
        sim.box_size = common.Vec(10, 10, 10)
        np.testing.assert_equal(sim.kernel_supports_topologies(), True)
        sim.register_particle_type("T", 1.0, .5, flavor=ParticleTypeFlavor.TOPOLOGY)
        p1 = sim.create_topology_particle("T", common.Vec(0, 0, 0))
        p2 = sim.create_topology_particle("T", common.Vec(1, 1, 1))
        p3 = sim.create_topology_particle("T", common.Vec(2, 1, 1))
        p4 = sim.create_topology_particle("T", common.Vec(3, 1, 1))
        np.testing.assert_equal(p1.get_position(), common.Vec(0, 0, 0))
        np.testing.assert_equal(p2.get_position(), common.Vec(1, 1, 1))
        np.testing.assert_equal(p3.get_position(), common.Vec(2, 1, 1))
        np.testing.assert_equal(p4.get_position(), common.Vec(3, 1, 1))
        particles = [p1, p2, p3, p4]
        top = sim.add_topology("TA", particles)
        bond1 = HarmonicBondPotentialBond(0, 1, 2.5, 3.3)
        bond2 = HarmonicBondPotentialBond(0, 2, 2.5, 3.3)
        top.add_harmonic_bond_potential([bond1, bond2])
        forces = []

        def callback(result):
            forces.append(result)

        sim.register_observable_forces(1, ["T"], callback)
        self.assertRaises(ValueError, sim.run_scheme_readdy(True).configure_and_run, 0, 1)

if __name__ == '__main__':
    unittest.main()
