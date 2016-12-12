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


from __future__ import print_function
import unittest

import readdy._internal.prototyping as pr
import readdy._internal.common as cmn

import numpy as np


class TestModel(unittest.TestCase):
    def setUp(self):
        self.kernel = pr.SingleCPUKernel()
        self.ctx = self.kernel.get_kernel_context()
        self.model = self.kernel.get_kernel_state_model()
        self.progs = self.kernel.get_program_factory()
        self.pots = self.kernel.get_potential_factory()

    def test_kernel_context_kbt(self):
        self.ctx.kbt = 5.0
        np.testing.assert_equal(self.ctx.kbt, 5.0)

    def test_kernel_context_box_size(self):
        self.ctx.set_box_size(cmn.Vec(5, 5, 5))
        np.testing.assert_equal(self.ctx.get_box_size(), cmn.Vec(5, 5, 5))

    def test_kernel_context_periodic_boundary(self):
        self.ctx.periodic_boundary = [True, False, True]
        np.testing.assert_equal(self.ctx.periodic_boundary, [True, False, True])

    def test_kernel_context_set_diffusion_constant(self):
        self.ctx.set_diffusion_constant("A", 13.0)
        np.testing.assert_equal(self.ctx.get_diffusion_constant("A"), 13.0)

    def test_kernel_context_fix_position_fun(self):
        self.ctx.set_box_size(cmn.Vec(1, 1, 1))
        self.ctx.periodic_boundary = [True, True, True]
        fix_pos = self.ctx.get_fix_position_fun()
        v_outside = cmn.Vec(4, 4, 4)
        fix_pos(v_outside)
        v_inside = cmn.Vec(.1, .1, .1)
        fix_pos(v_inside)
        np.testing.assert_equal(v_outside, cmn.Vec(0, 0, 0))
        np.testing.assert_equal(v_inside, cmn.Vec(.1, .1, .1))

    def test_kernel_context_shortest_difference(self):
        self.ctx.set_box_size(cmn.Vec(2, 2, 2))
        self.ctx.periodic_boundary = [True, True, True]
        diff = self.ctx.get_shortest_difference_fun()
        np.testing.assert_equal(diff(cmn.Vec(0, 0, 0), cmn.Vec(1, 0, 0)), cmn.Vec(1, 0, 0),
                          err_msg="both vectors were already inside the domain")
        np.testing.assert_equal(diff(cmn.Vec(0, 0, 0), cmn.Vec(-1.5, 0, 0)), cmn.Vec(.5, 0, 0),
                          err_msg="The latter vector was outside the domain, thus had to get a "
                                  "position update to (.5, 0, 0)")

    def test_kernel_context_dist_squared_fun(self):
        self.ctx.set_box_size(cmn.Vec(2, 2, 2))
        self.ctx.periodic_boundary = [True, True, True]
        dist = self.ctx.get_dist_squared_fun()
        np.testing.assert_equal(dist(cmn.Vec(0, 0, 0), cmn.Vec(1, 0, 0)), 1.0)
        np.testing.assert_equal(dist(cmn.Vec(-.5, 0, 0), cmn.Vec(-1.5, 0, 0)), 1.0)

    def test_kernel_context_particle_radius(self):
        self.ctx.set_diffusion_constant("A", 1.0)
        self.ctx.set_particle_radius("A", 5.0)
        np.testing.assert_equal(self.ctx.get_particle_radius("A"), 5.0)

    def test_kernel_context_time_step(self):
        self.ctx.timestep = 111.0
        np.testing.assert_equal(self.ctx.timestep, 111.0)

    def test_potential_order_1(self):
        self.ctx.set_box_size(cmn.Vec(2, 2, 2))
        self.ctx.periodic_boundary = [True, True, True]
        self.model.get_particle_data().clear()
        class MyPot1(pr.PotentialOrder1):
            def __init__(self):
                super(MyPot1, self).__init__("yes")

            def calculate_energy(self, pos_vec):
                return 5.0

            def calculate_force(self, pos_vec):
                return cmn.Vec(.1, .1, .1)

            def configure_for_type(self, type):
                pass

            def get_relevant_length_scale(self):
                return 5.0

            def get_maximal_force(self, kbt):
                return kbt

        self.ctx.set_diffusion_constant("A", 1.0)
        pot = MyPot1()
        self.ctx.register_potential_order_1(pot, "A")
        add_particles_program = self.progs.create_add_particles()
        add_particles_program.add_particle(pr.Particle(0, 0, .5, self.ctx.get_particle_type_id("A"))) # x y z type_id
        add_particles_program.execute()
        self.ctx.configure()
        updforces = self.progs.create_update_forces()
        updforces.execute()

        np.testing.assert_equal(self.model.get_energy(), 5.0, err_msg="the user defined potential returns energy=5.0")

        data = self.model.get_particle_data()
        it_forces = data.forces
        f_vec = next(it_forces)
        np.testing.assert_equal(f_vec, cmn.Vec(.1, .1, .1))
        with np.testing.assert_raises(StopIteration):
            next(it_forces)

    def test_potential_order_2(self):
        self.ctx.set_box_size(cmn.Vec(2, 2, 2))
        self.ctx.periodic_boundary = [True, True, True]
        self.model.get_particle_data().clear()

        class MyPot2(pr.PotentialOrder2):
            def __init__(self):
                super(MyPot2, self).__init__("yes2")
            
            def calculate_energy(self, x_ij):
                return np.sqrt(x_ij * x_ij)

            def calculate_force(self, x_ij):
                return .5 * x_ij

            def configure_for_types(self, type1, type2):
                pass

            def get_cutoff_radius(self):
                return 5.0

            def get_maximal_force(self, kbt):
                return kbt

        self.ctx.set_diffusion_constant("A", 1.0)
        self.ctx.set_diffusion_constant("B", 1.0)
        pot = MyPot2()
        self.ctx.register_potential_order_2(pot, "A", "B")
        add_particles_program = self.progs.create_add_particles()
        add_particles_program.add_particle(pr.Particle(0, 0, 0, self.ctx.get_particle_type_id("A")))
        add_particles_program.add_particle(pr.Particle(1, 1, 1, self.ctx.get_particle_type_id("B")))
        add_particles_program.execute()
        self.ctx.configure()
        self.progs.create_update_neighbor_list().execute()
        self.progs.create_update_forces().execute()

        np.testing.assert_almost_equal(self.model.get_energy(), np.sqrt(3))

        it_forces = self.model.get_particle_data().forces
        it_types = self.model.get_particle_data().types
        t = next(it_types)
        if t == self.ctx.get_particle_type_id("A"):
            np.testing.assert_equal(next(it_forces), cmn.Vec(.5, .5, .5))
            np.testing.assert_equal(next(it_forces), cmn.Vec(-.5, -.5, -.5))
        else:
            np.testing.assert_equal(next(it_forces), cmn.Vec(-.5, -.5, -.5))
            np.testing.assert_equal(next(it_forces), cmn.Vec(.5, .5, .5))

        with np.testing.assert_raises(StopIteration):
            next(it_forces)

    def test_potential_factory(self):
        cube_pot = self.pots.create_cube_potential()
        np.testing.assert_equal(cube_pot.get_name(), "Cube")

        repulsion_pot = self.pots.create_harmonic_repulsion()
        np.testing.assert_equal(repulsion_pot.get_name(), "HarmonicRepulsion")

        weak_interaction_pot = self.pots.create_weak_interaction()
        np.testing.assert_equal(weak_interaction_pot.get_name(), "WeakInteractionPiecewiseHarmonic")

if __name__ == '__main__':
    unittest.main()
