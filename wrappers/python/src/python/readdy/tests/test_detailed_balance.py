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
Created on 04.12.17

@author: chrisfroe
"""

import readdy
from readdy.util.testing_utils import ReaDDyTestCase


class TestDetailedBalanceReactionHandler(ReaDDyTestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()

    @classmethod
    def tearDownClass(cls):
        super().tearDownClass()

    def test_low_concentration_limit(self):
        rds = readdy.ReactionDiffusionSystem(box_size=[10., 10., 10.], unit_system=None)
        rds.add_species("A", 1.)
        rds.add_species("B", 1.)
        rds.add_species("C", 1.)
        rds.reactions.add("fus: A +(3) B -> C", 10.)
        rds.reactions.add("fis: C -> A +(3) B", 2.)
        rds.reactions.add("conversion: A -> C", 5.)
        rds.potentials.add_harmonic_repulsion("A", "B", 2., 3.)
        simulation = rds.simulation("SingleCPU")
        simulation.reaction_handler = "DetailedBalance"
        simulation.add_particle("A", [0.,0.,0.])
        simulation.add_particle("B", [0.,0.,0.])
        simulation.run(1000, 1e-10, show_system=False)
