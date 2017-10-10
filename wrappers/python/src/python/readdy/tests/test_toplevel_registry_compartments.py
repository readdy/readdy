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
@author: chrisfroe
"""

from readdy.api.reaction_diffusion_system import ReactionDiffusionSystem
from readdy.util.testing_utils import ReaDDyTestCase


class TestToplevelRegistryCompartments(ReaDDyTestCase):
    def setUp(self):
        self.rds = ReactionDiffusionSystem(box_size=[1., 1., 1.])
        self.rds.add_species("A", 1.)
        self.rds.add_species("B", 1.)

    def test_sphere_sanity(self):
        self.rds.compartments.add_sphere({"A": "B"}, "mysphere", [0, 1, 2], 2., False)

    def test_sphere_invalid_args(self):
        with self.assertRaises(ValueError):
            self.rds.compartments.add_sphere({"A": "B"}, "mysphere", [0, 1, 2, 4], 2., False)

        with self.assertRaises(ValueError):
            self.rds.compartments.add_sphere(["A", "B"], "mysphere", [0, 1, 2], 2., False)

        with self.assertRaises(ValueError):
            self.rds.compartments.add_sphere({"A": "B"}, "mysphere", [0, 1, 2], 0., False)

        with self.assertRaises(ValueError):
            self.rds.compartments.add_sphere({"A": "B"}, "mysphere", [0, 1, 2], 1., "nobool")

    def test_plane_sanity(self):
        self.rds.compartments.add_plane({"A": "B"}, "myplane", [0, 1, 0], 2., False)

    def test_plane_invalid_args(self):
        with self.assertRaises(ValueError):
            self.rds.compartments.add_plane({"A": "B"}, "myplane", [0, 1, 0, 0], 2., False)

        with self.assertRaises(ValueError):
            self.rds.compartments.add_plane(["A", "B"], "myplane", [0, 1, 0], 2., False)

        with self.assertRaises(ValueError):
            self.rds.compartments.add_plane({"A": "B"}, "myplane", [0, 1, 0], -1., False)

        with self.assertRaises(ValueError):
            self.rds.compartments.add_plane({"A": "B"}, "myplane", [0, 1, 0], 2., "nobool")
