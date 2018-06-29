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
from readdy._internal.readdybinding.api import Simulation

from readdy.util.testing_utils import ReaDDyTestCase


class TestTimer(ReaDDyTestCase):

    def test_timer_sanity(self):
        simulation = Simulation("CPU")
        simulation.run(10, .1)
        root = simulation.performance_root()
        np.testing.assert_equal(root.count(), 1)
        np.testing.assert_equal(root.time() > 0., True)

if __name__ == '__main__':
    unittest.main()
