# coding=utf-8

# Copyright © 2017 Computational Molecular Biology Group,
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

import unittest
import numpy as np
import readdy._internal.readdybinding.common.perf as perf
from readdy._internal.readdybinding.api import Simulation

class TestTimer(unittest.TestCase):
    def test_timer_sanity(self):
        perf.clear()
        simulation = Simulation()
        simulation.set_kernel("SingleCPU")
        scheme = simulation.run_scheme_readdy(True)
        scheme.configure_and_run(10, 0.1)
        np.testing.assert_equal(perf.counts()["integrator"], 10)
        np.testing.assert_equal(perf.times()["integrator"] >= 0., True)


if __name__ == '__main__':
    unittest.main()