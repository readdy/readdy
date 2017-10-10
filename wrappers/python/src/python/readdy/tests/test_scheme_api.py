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
import numpy as np

from readdy._internal.readdybinding.api import Simulation
from readdy._internal.readdybinding.common import Vec

from readdy.util.testing_utils import ReaDDyTestCase


class TestSchemeApi(ReaDDyTestCase):
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
            .configure_and_run(10, 1)

        simulation.run_scheme_readdy(False) \
            .with_integrator("EulerBDIntegrator") \
            .include_forces(False) \
            .with_reaction_scheduler("UncontrolledApproximation") \
            .evaluate_observables(False) \
            .configure(1).run(10)

    def test_interrupt_simple(self):
        sim = Simulation()
        sim.set_kernel("SingleCPU")
        sim.register_particle_type("A", 0.1)
        # Define counter as list. This is a workaround because nosetest will complain otherwise.
        counter = [0]

        def increment(result):
            counter[0] += 1

        sim.register_observable_n_particles(1, ["A"], increment)
        scheme = sim.run_scheme_readdy(True).configure(0.1)
        do_continue = lambda t: t < 5
        scheme.run_with_criterion(do_continue)
        np.testing.assert_equal(counter[0], 6)

    def test_interrupt_maxparticles(self):
        sim = Simulation()
        sim.set_kernel("SingleCPU")
        sim.register_particle_type("A", 0.1)
        sim.add_particle("A", Vec(0, 0, 0))
        sim.register_reaction_fission("bla", "A", "A", "A", 1000., 0., 0.5, 0.5)
        counter = [0]
        shall_stop = [False]

        def increment(result):
            counter[0] += 1
            if result[0] >= 8:
                shall_stop[0] = True

        sim.register_observable_n_particles(1, ["A"], increment)
        scheme = sim.run_scheme_readdy(True).configure(1.)
        do_continue = lambda t: not shall_stop[0]
        scheme.run_with_criterion(do_continue)
        np.testing.assert_equal(counter[0], 4)


if __name__ == '__main__':
    unittest.main()
