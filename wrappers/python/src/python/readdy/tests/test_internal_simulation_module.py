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
import unittest

import numpy as np
from readdy._internal.readdybinding.api import Pot2
from readdy._internal.readdybinding.api import Simulation
from readdy._internal.readdybinding.api import KernelProvider

from readdy._internal.readdybinding.common import Vec

from readdy.util import platform_utils
from readdy.util.testing_utils import ReaDDyTestCase


class TestInternalSimulationModule(ReaDDyTestCase):

    def setUp(self):
        super().setUp()
        self.kernel_provider = KernelProvider.get()
        self.kernel_provider.load_from_dir(platform_utils.get_readdy_plugin_dir())
        self.simulation = Simulation()

    def py_harmonic_repulsion_energy(self, x_ij):
        dist = x_ij * x_ij
        # if dist < sqrt(25): return energy with  force constant 1
        if dist < 25:
            return (np.sqrt(dist) - 5) ** 2
        else:
            return 0

    def py_harmonic_repulsion_force(self, x_ij):
        dist = x_ij * x_ij
        if dist < 25:
            dist = np.sqrt(dist)
            return (2 * (dist - 5) / dist) * x_ij
        else:
            return Vec(0, 0, 0)

    def test_properties(self):
        if not self.simulation.is_kernel_selected():
            self.simulation.set_kernel('SingleCPU')
        np.testing.assert_equal(self.simulation.is_kernel_selected(), True)
        np.testing.assert_equal(self.simulation.get_selected_kernel_type(), "SingleCPU")
        self.simulation.kbt = 5.0
        np.testing.assert_equal(self.simulation.kbt, 5.0)
        self.simulation.periodic_boundary = [True, False, True]
        np.testing.assert_equal(self.simulation.periodic_boundary, (True, False, True))
        self.simulation.box_size = Vec(1, 3.6, 7)
        np.testing.assert_equal(self.simulation.box_size, Vec(1, 3.6, 7))

    def py_position_observable_callback(self, positions):
        _vec_sum = Vec(0, 0, 0)
        for v in positions:
            _vec_sum += v
        mean = _vec_sum / float(len(positions))
        print("mean=%s" % mean)

    def test_potentials(self):
        if not self.simulation.is_kernel_selected():
            self.simulation.set_kernel("SingleCPU")

        ida = self.simulation.register_particle_type("ParticleTypeA", 1.0)
        idb = self.simulation.register_particle_type("ParticleTypeB", 3.0)
        self.simulation.register_particle_type("ParticleTypeA_internal", 1.0)
        self.simulation.register_particle_type("ParticleTypeB_internal", 3.0)
        pot = Pot2(ida, idb, self.py_harmonic_repulsion_energy,
                   self.py_harmonic_repulsion_force)
        self.simulation.register_potential_order_2(pot)
        self.simulation.register_potential_harmonic_repulsion("ParticleTypeA_internal", "ParticleTypeB_internal", 1.0, 2.0)
        self.simulation.add_particle("ParticleTypeA", Vec(0, 0, 0))
        self.simulation.add_particle("ParticleTypeB", Vec(0.4, 0.4, 0.4))
        self.simulation.run(100, 1)


if __name__ == '__main__':
    unittest.main()
