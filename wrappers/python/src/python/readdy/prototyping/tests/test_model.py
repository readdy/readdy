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

    def test_kernel_context_register_particle_type(self):
        self.ctx.register_particle_type("A", 13.0, 17.0)
        np.testing.assert_equal(self.ctx.get_diffusion_constant("A"), 13.0)
        np.testing.assert_equal(self.ctx.get_particle_radius("A"), 17.0)

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

    def test_potential_order_1(self):
        self.ctx.set_box_size(cmn.Vec(2, 2, 2))
        self.ctx.periodic_boundary = [True, True, True]
        self.model.get_particle_data().clear()
        class MyPot1(pr.PotentialOrder1):
            def __init__(self, type):
                super(MyPot1, self).__init__(type)

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

        self.ctx.register_particle_type("A", 1.0, 1.0)
        pot = MyPot1("A")
        self.ctx.register_potential_order_1(pot)
        particles = [pr.Particle(0, 0, .5, self.ctx.get_particle_type_id("A"))]
        add_particles_program = self.progs.create_add_particles(particles)
        add_particles_program.perform()
        self.ctx.configure()
        updforces = self.progs.create_update_forces()
        updforces.perform()

        np.testing.assert_equal(self.model.get_energy(), 5.0, err_msg="the user defined potential returns energy=5.0")

        data = self.model.get_particle_data()
        it = data.entries
        f_vec = next(it).force
        np.testing.assert_equal(f_vec, cmn.Vec(.1, .1, .1))
        with np.testing.assert_raises(StopIteration):
            next(it)

    def test_potential_order_2(self):
        self.ctx.set_box_size(cmn.Vec(2, 2, 2))
        self.ctx.periodic_boundary = [True, True, True]
        self.model.get_particle_data().clear()

        class MyPot2(pr.PotentialOrder2):
            def __init__(self, typea, typeb):
                super(MyPot2, self).__init__(typea, typeb)
            
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

        self.ctx.register_particle_type("A", 1.0, 1.0)
        self.ctx.register_particle_type("B", 1.0, 1.0)
        pot = MyPot2("A", "B")
        self.ctx.register_potential_order_2(pot)
        particles = [pr.Particle(0, 0, 0, self.ctx.get_particle_type_id("A")),
                     pr.Particle(1, 1, 1, self.ctx.get_particle_type_id("B"))]
        add_particles_program = self.progs.create_add_particles(particles)
        add_particles_program.perform()
        self.ctx.configure()
        self.progs.create_update_neighbor_list().perform()
        self.progs.create_update_forces().perform()

        np.testing.assert_almost_equal(self.model.get_energy(), np.sqrt(3))

        it = self.model.get_particle_data().entries
        entry = next(it)
        if entry.type == self.ctx.get_particle_type_id("A"):
            np.testing.assert_equal(entry.force, cmn.Vec(.5, .5, .5))
            np.testing.assert_equal(next(it).force, cmn.Vec(-.5, -.5, -.5))
        else:
            np.testing.assert_equal(entry.force, cmn.Vec(-.5, -.5, -.5))
            np.testing.assert_equal(next(it).force, cmn.Vec(.5, .5, .5))

        with np.testing.assert_raises(StopIteration):
            next(it)

    def test_potential_factory(self):
        self.ctx.register_particle_type("A", 1.0, 1.0)
        self.ctx.register_particle_type("B", 1.0, 1.0)
        cube_pot = self.pots.create_cube_potential("A", 1, cmn.Vec(0,0,0), cmn.Vec(1,1,1), True)
        np.testing.assert_equal(cube_pot.get_name(), "Cube")

        repulsion_pot = self.pots.create_harmonic_repulsion("A", "B", 10)
        np.testing.assert_equal(repulsion_pot.get_name(), "HarmonicRepulsion")

        weak_interaction_pot = self.pots.create_weak_interaction("A", "B", 0,0,0,0)
        np.testing.assert_equal(weak_interaction_pot.get_name(), "WeakInteractionPiecewiseHarmonic")

if __name__ == '__main__':
    unittest.main()
