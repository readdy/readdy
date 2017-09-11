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
Created on 08.09.17

@author: clonker
"""
from readdy.api.observables import Observables


class Simulation(object):

    def __init__(self, reaction_diffusion_system, kernel):
        self._reaction_diffusion_system = reaction_diffusion_system
        self._kernel = kernel
        self._output_file = ""
        self._observables = Observables(self)
        self._integrator = "EulerBDIntegrator"
        self._reaction_handler = "Gillespie"
        # todo this should probably (default) determined depending on whether top reactions were defined
        self._evaluate_topology_reactions = True
        

    @property
    def output_file(self):
        return self._output_file

    @output_file.setter
    def output_file(self, value):
        self._output_file = value

    @property
    def observe(self):
        return self._observables