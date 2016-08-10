import unittest

import readdy._internal.prototyping as pr
import readdy._internal.common as cmn
import readdy._internal.api as api

import numpy.testing as test


class TestModel(unittest.TestCase):
    def setUp(self):
        self.kernel = pr.SingleCPUKernel()
        self.ctx = self.kernel.get_kernel_context()
        self.model = self.kernel.get_kernel_state_model()
        self.progs = self.kernel.get_program_factory()

    def test_kernel_context_kbt(self):
        self.ctx.kbt = 5.0
        test.assert_equal(self.ctx.kbt, 5.0)

    def test_kernel_context_box_size(self):
        self.ctx.box_size = cmn.Vec(5, 5, 5)
        test.assert_equal(self.ctx.box_size, cmn.Vec(5, 5, 5))

    def test_kernel_context_periodic_boundary(self):
        self.ctx.periodic_boundary = [True, False, True]
        test.assert_equal(self.ctx.periodic_boundary, [True, False, True])

    def test_kernel_context_set_diffusion_constant(self):
        self.ctx.set_diffusion_constant("A", 13.0)
        test.assert_equal(self.ctx.get_diffusion_constant("A"), 13.0)

    def test_kernel_context_fix_position_fun(self):
        self.ctx.box_size = cmn.Vec(1, 1, 1)
        self.ctx.periodic_boundary = [True, True, True]
        fix_pos = self.ctx.get_fix_position_fun()
        v_outside = cmn.Vec(4, 4, 4)
        fix_pos(v_outside)
        v_inside = cmn.Vec(.1, .1, .1)
        fix_pos(v_inside)
        test.assert_equal(v_outside, cmn.Vec(0, 0, 0))
        test.assert_equal(v_inside, cmn.Vec(.1, .1, .1))

    def test_kernel_context_shortest_difference(self):
        self.ctx.box_size = cmn.Vec(2, 2, 2)
        self.ctx.periodic_boundary = [True, True, True]
        diff = self.ctx.get_shortest_difference_fun()
        test.assert_equal(diff(cmn.Vec(0, 0, 0), cmn.Vec(1, 0, 0)), cmn.Vec(1, 0, 0),
                          err_msg="both vectors were already inside the domain")
        test.assert_equal(diff(cmn.Vec(0, 0, 0), cmn.Vec(-1.5, 0, 0)), cmn.Vec(.5, 0, 0),
                          err_msg="The latter vector was outside the domain, thus had to get a "
                                  "position update to (.5, 0, 0)")

    def test_kernel_context_dist_squared_fun(self):
        self.ctx.box_size = cmn.Vec(2, 2, 2)
        self.ctx.periodic_boundary = [True, True, True]
        dist = self.ctx.get_dist_squared_fun()
        test.assert_equal(dist(cmn.Vec(0, 0, 0), cmn.Vec(1, 0, 0)), 1.0)
        test.assert_equal(dist(cmn.Vec(-.5, 0, 0), cmn.Vec(-1.5, 0, 0)), 1.0)

    def test_kernel_context_particle_radius(self):
        self.ctx.set_diffusion_constant("A", 1.0)
        self.ctx.set_particle_radius("A", 5.0)
        test.assert_equal(self.ctx.get_particle_radius("A"), 5.0)

    def test_kernel_context_time_step(self):
        self.ctx.timestep = 111.0
        test.assert_equal(self.ctx.timestep, 111.0)

    def test_potential_order_1(self):
        class MyPot1(pr.PotentialOrder1):
            def __init__(self):
                super(MyPot1, self).__init__("yes")

            def calculate_energy(self, pos_vec):
                return 5.0

            def calculate_force(self, f_vec, pos_vec):
                f_vec += cmn.Vec(.1, .1, .1)

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
        add_particles_program.add_particle(pr.Particle(0,0,.5,1))
        add_particles_program.execute()
        self.ctx.configure()
        update_forces = self.progs.create_update_forces()
        update_forces.execute()
        # todo: get particle data object and check if forces and energy are correct (according to this potential)


if __name__ == '__main__':
    unittest.main()
