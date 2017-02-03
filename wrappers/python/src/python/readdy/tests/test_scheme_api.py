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


import unittest

from readdy._internal.api import Simulation


class TestSchemeApi(unittest.TestCase):
    def test_sanity(self):
        simulation = Simulation()
        simulation.set_kernel("SingleCPU")
        configurator = simulation.run_scheme_readdy(False)
        scheme = configurator \
            .with_integrator("EulerBDIntegrator") \
            .include_forces(False) \
            .with_reaction_scheduler("UncontrolledApproximation") \
            .evaluate_observables(False) \
            .configure(1)
        scheme.run(10)

    def test_sanity_advanced(self):
        # todo @chrisfroe test interrupt api
        simulation = Simulation()
        simulation.set_kernel("SingleCPU")
        configurator = simulation.run_scheme_advanced(False)
        scheme = configurator \
            .with_integrator("EulerBDIntegrator") \
            .include_forces(False) \
            .include_compartments(False) \
            .with_reaction_scheduler("UncontrolledApproximation") \
            .evaluate_observables(False) \
            .configure(1)
        scheme.run(10)

    def test_sanity_oneliner(self):
        simulation = Simulation()
        simulation.set_kernel("SingleCPU")

        simulation.run_scheme_readdy(False) \
            .with_integrator("EulerBDIntegrator") \
            .include_forces(False) \
            .with_reaction_scheduler("UncontrolledApproximation") \
            .evaluate_observables(False) \
            .configure_and_run(1, 10)

        simulation.run_scheme_readdy(False) \
            .with_integrator("EulerBDIntegrator") \
            .include_forces(False) \
            .with_reaction_scheduler("UncontrolledApproximation") \
            .evaluate_observables(False) \
            .configure(1).run(10)


if __name__ == '__main__':
    unittest.main()
