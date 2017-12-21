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
Created on 10.10.17

@author: clonker
"""

import readdy
from readdy.util.testing_utils import ReaDDyTestCase

from readdy.api.conf.UnitConfiguration import UnitConfiguration

import numpy as np
import pint

class TestUnits(ReaDDyTestCase):

    def test_unit_configuration_sanity(self):
        conf = UnitConfiguration()
        reg = pint.UnitRegistry()
        np.testing.assert_equal(1 * conf.length_unit, 1 * reg.nanometer)
        np.testing.assert_equal(1 * conf.time_unit, 1 * reg.nanoseconds)
        np.testing.assert_equal(1 * conf.energy_unit, 1 * reg.kilojoule / reg.mol)

    def test_custom_configuration(self):
        conf = UnitConfiguration(length_unit='nanometer', time_unit='millisecond', energy_unit='kilocal/mol')
        a = 1*conf.length_unit
        b = (1e-6*conf.reg.millimeter).to(conf.length_unit)
        np.testing.assert_almost_equal(a.magnitude, b.magnitude)

    def test_force_constant_unit(self):
        conf = UnitConfiguration()
        kunit = conf.energy_unit / conf.length_unit**2
        np.testing.assert_equal(kunit, conf.reg.parse_units("kilojoule/(mol * nanometer * nanometer)"))

    def test_box_size(self):
        sys = readdy.ReactionDiffusionSystem([5,5,5] * readdy.units.meter)
        np.testing.assert_equal(sys.box_size.magnitude, ([5, 5, 5] * sys.units.meter).to(sys.length_unit).magnitude)

    def test_radians(self):
        sys = readdy.ReactionDiffusionSystem(box_size=[1., 1., 1.])
        angle = 180 * sys.units.degree
        np.testing.assert_equal(np.pi, angle.to(sys.units.radians).magnitude)
