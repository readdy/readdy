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
