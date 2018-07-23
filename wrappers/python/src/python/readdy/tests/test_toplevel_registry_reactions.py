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
Created on 27.09.17

@author: chrisfroe
"""

from readdy.api.reaction_diffusion_system import ReactionDiffusionSystem
from readdy.util.testing_utils import ReaDDyTestCase

import readdy as _readdy
ut = _readdy.units


class TestToplevelRegistryReactions(ReaDDyTestCase):
    def setUp(self):
        self.rds = ReactionDiffusionSystem(box_size=[1., 1., 1.])
        self.rds.add_species("A", 1.)
        self.rds.add_species("B", 1.)
        self.rds.add_species("C", 1.)

    def test_add_decay(self):
        self.rds.reactions.add_decay("label1", "A", 1.0)
        self.rds.reactions.add_decay("label2", "A", 1.0 / ut.second)

    def test_add_conversion(self):
        self.rds.reactions.add_conversion("label1", "A", "B", 2.0)
        self.rds.reactions.add_conversion("label2", "A", "B", 2.0 / ut.second)

    def test_add_fusion(self):
        self.rds.reactions.add_fusion("label1", "A", "B", "C", 2.0, 3.0)
        self.rds.reactions.add_fusion("label2", "A", "B", "C", 2.0 * 1. / ut.second, 3.0 * ut.nanometer, weight1=0., weight2=1.)

    def test_add_fission(self):
        self.rds.reactions.add_fission("label", "C", "A", "B", 2.0, 3.0, weight1=0., weight2=1.)
        self.rds.reactions.add_fission("label2", "C", "A", "B", 2.0 * 1. / ut.second, 3.0 * ut.nanometer)

    def test_add_enzymatic(self):
        self.rds.reactions.add_enzymatic("label1", "C", "A", "B", 2.0, 4.0)
        self.rds.reactions.add_enzymatic("label2", "C", "A", "B", 2.0 / ut.second, 4.0 * ut.nanometer)

    def test_add_invalid_descriptor(self):
        with self.assertRaises(ValueError):
            self.rds.reactions.add("enz: A +(1) B -> A + A", 5.)

    def test_add_descriptor(self):
        self.rds.reactions.add("label1: A +(1) B -> C", 5.)
        self.rds.reactions.add("label2: A +(1) B -> C", 5. / ut.nanosecond)
