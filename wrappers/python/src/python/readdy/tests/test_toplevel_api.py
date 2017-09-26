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
Created on 26.09.17

@author: clonker
"""

import numpy as np

from readdy.api.reaction_diffusion_system import ReactionDiffusionSystem
from readdy.util.testing_utils import ReaDDyTestCase


class TestTopologies(ReaDDyTestCase):

    def test_kbt(self):
        rdf = ReactionDiffusionSystem()
        rdf.kbt = 5.
        np.testing.assert_equal(rdf.kbt, 5.)

    def test_box_size(self):
        rdf = ReactionDiffusionSystem()
        rdf.box_size = [1., 2., 3.]
        np.testing.assert_equal(rdf.box_size, [1., 2., 3.])
        rdf.box_size = np.array([5., 6., 7.])
        np.testing.assert_equal(rdf.box_size, [5., 6., 7.])
        rdf.box_size = (1., 5., 7.)
        np.testing.assert_equal(rdf.box_size, [1., 5., 7.])
        np.testing.assert_equal(rdf.box_volume, 5.*7.)

    def test_pbc(self):
        rdf = ReactionDiffusionSystem()
        rdf.periodic_boundary_conditions = True, False, True
        np.testing.assert_equal(rdf.periodic_boundary_conditions, [True, False, True])
        rdf.periodic_boundary_conditions = np.array([False, False, True])
        np.testing.assert_equal(rdf.periodic_boundary_conditions, [False, False, True])
