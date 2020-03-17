# coding=utf-8

# Copyright © 2020 Computational Molecular Biology Group,
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
Created on 17.03.20

@author: chrisfroe
"""

import os
import readdy
import unittest
import numpy as np
import tempfile
import shutil
from tqdm import tqdm as tqdm


class TestCustomLoop(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.dir = tempfile.mkdtemp("test-custom-loop")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.dir, ignore_errors=True)

    def test_displace_particle(self):
        rds = readdy.ReactionDiffusionSystem(box_size=[10., 10., 10.], unit_system=None)
        rds.add_species("A", 1.0)
        simulation = rds.simulation("SingleCPU")
        simulation.add_particle("A", [0., 0., 0.])
        dt = 1.0
        # action factory is an experimental feature and thus _hidden
        integrator = simulation._actions.integrator_euler_brownian_dynamics(dt)
        p = simulation.current_particles[0]
        self.assertTrue(np.all(p.pos == np.array([0., 0., 0.])))
        integrator()
        p = simulation.current_particles[0]
        self.assertFalse(np.all(p.pos == np.array([0., 0., 0.])))

    def test_build_and_run_custom_loop(self):
        rds = readdy.ReactionDiffusionSystem(box_size=[7., 7., 7.], unit_system=None)
        rds.add_species("A", 1.0)
        rds.add_species("B", 1.0)
        rds.add_species("C", 1.0)
        rds.reactions.add("fusion: A +(2) B -> C", 10.)
        simulation = rds.simulation("SingleCPU")
        simulation.add_particle("A", [0., 0., 0.])
        simulation.add_particle("B", [3., 0., 0.])
        simulation.output_file = os.path.join(self.dir, "customlooptest1.h5")
        simulation.observe.number_of_particles(1, ["A", "B", "C"])

        def loop():
            nonlocal simulation
            dt = 1.0
            n_steps = 10000
            base_path = os.path.join(self.dir, "ckpts")
            os.makedirs(base_path, exist_ok=True)
            max_n_saves = 2

            diff = simulation._actions.integrator_euler_brownian_dynamics(dt)
            calc_forces = simulation._actions.calculate_forces()
            create_nl = simulation._actions.create_neighbor_list(rds.calculate_max_cutoff())
            update_nl = simulation._actions.update_neighbor_list()
            reac = simulation._actions.reaction_handler_uncontrolled_approximation(dt)
            obs = simulation._actions.evaluate_observables()
            check = simulation._actions.make_checkpoint(base_path, max_n_saves)

            create_nl()
            calc_forces()
            update_nl()
            obs(0)
            for t in tqdm(range(1, n_steps + 1)):
                diff()
                update_nl()
                reac()
                update_nl()
                calc_forces()
                obs(t)  # striding of observables is done internally
                if t % 100 == 0:
                    check(t)

        simulation.run_custom_loop(loop)

        traj = readdy.Trajectory(simulation.output_file)
        ts, ns = traj.read_observable_number_of_particles()

        self.assertEqual(ns[0, 0], 1)
        self.assertEqual(ns[0, 1], 1)
        self.assertEqual(ns[0, 2], 0)
        self.assertEqual(ns[-1, 0], 0)
        self.assertEqual(ns[-1, 1], 0)
        self.assertEqual(ns[-1, 2], 1)
