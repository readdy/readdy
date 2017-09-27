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
from readdy.api.conf.KernelConfiguration import CPUKernelConfiguration as _CPUKernelConfiguration
from readdy.api.conf.KernelConfiguration import NOOPKernelConfiguration as _NOOPKernelConfiguration
from readdy.api.registry.observables import Observables as _Observables
from readdy._internal.readdybinding.api import Simulation as _Simulation

class Simulation(object):

    def __init__(self, kernel, context):
        self._kernel = kernel
        self._simulation = _Simulation()
        self._simulation.set_kernel(kernel)
        self._simulation.context = context

        self._output_file = ""
        self._observables = _Observables(self)
        self._integrator = "EulerBDIntegrator"
        self._reaction_handler = "Gillespie"
        # todo this should probably (default) determined depending on whether top reactions were defined
        self._evaluate_topology_reactions = True
        self._simulation_scheme = "DefaultScheme"

        if kernel == "CPU":
            self._kernel_configuration = _CPUKernelConfiguration()
        else:
            self._kernel_configuration = _NOOPKernelConfiguration()
        

    @property
    def simulation_scheme(self):
        return self._simulation_scheme

    @simulation_scheme.setter
    def simulation_scheme(self, value):
        if value in ("DefaultScheme", "AdvancedScheme"):
            self._simulation_scheme = value
        else:
            raise ValueError("Simulation scheme value can only be one of \"DefaultScheme\" and \"AdvancedScheme\".")

    @property
    def kernel_configuration(self):
        return self._kernel_configuration

    @property
    def output_file(self):
        return self._output_file

    @output_file.setter
    def output_file(self, value):
        self._output_file = value

    @property
    def observe(self):
        return self._observables