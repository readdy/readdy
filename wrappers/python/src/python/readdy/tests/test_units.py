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
