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


import unittest
import numpy as np

from readdy._internal.readdybinding.api import Simulation
from readdy._internal.readdybinding.common import Vec

from readdy.util.testing_utils import ReaDDyTestCase


class TestSchemeApi(ReaDDyTestCase):
    def test_sanity(self):
        simulation = Simulation("SingleCPU")
        loop = simulation.create_loop(1.)
        loop.use_integrator("EulerBDIntegrator")
        loop.evaluate_forces(False)
        loop.use_reaction_scheduler("UncontrolledApproximation")
        loop.evaluate_observables(False)
        loop.run(10)

    def test_interrupt_simple(self):
        sim = Simulation("SingleCPU")
        sim.context.particle_types.add("A", 0.1)
        # Define counter as list. This is a workaround because nosetest will complain otherwise.
        counter = [0]

        def increment(result):
            counter[0] += 1

        sim.register_observable_n_particles(1, ["A"], increment)
        do_continue = lambda t: t < 5
        sim.create_loop(.1).run_with_criterion(do_continue)
        np.testing.assert_equal(counter[0], 6)

    def test_interrupt_maxparticles(self):
        sim = Simulation("SingleCPU")
        sim.context.particle_types.add("A", 0.1)
        sim.add_particle("A", Vec(0, 0, 0))
        sim.context.reactions.add_fission("bla", "A", "A", "A", 1000., 0., 0.5, 0.5)
        counter = [0]
        shall_stop = [False]

        def increment(result):
            counter[0] += 1
            if result[0] >= 8:
                shall_stop[0] = True

        sim.register_observable_n_particles(1, ["A"], increment)
        loop = sim.create_loop(1.)
        do_continue = lambda t: not shall_stop[0]
        loop.run_with_criterion(do_continue)
        np.testing.assert_equal(counter[0], 4)


if __name__ == '__main__':
    unittest.main()
