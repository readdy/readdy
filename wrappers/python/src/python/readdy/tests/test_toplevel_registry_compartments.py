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
