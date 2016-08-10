import unittest

import readdy._internal.prototyping as pr
import readdy._internal.common as cmn
import readdy._internal.api as api

import numpy.testing as test


class TestModel(unittest.TestCase):
    def setUp(self):
        self.kernel = pr.SingleCPUKernel()
        self.ctx = self.kernel.get_kernel_context()

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
        self.ctx.box_size = cmn.Vec(1,1,1)
        self.ctx.periodic_boundary = [True, True, True]
        fix_pos = self.ctx.get_fix_position_fun()
        v_outside = cmn.Vec(4,4,4)
        fix_pos(v_outside)
        v_inside = cmn.Vec(.1, .1, .1)
        fix_pos(v_inside)
        test.assert_equal(v_outside, cmn.Vec(0,0,0))
        test.assert_equal(v_inside, cmn.Vec(.1, .1, .1))

    def test_kernel_context_shortest_difference(self):
        self.ctx.box_size = cmn.Vec(2,2,2)
        self.ctx.periodic_boundary = [True, True, True]
        diff = self.ctx.get_shortest_difference_fun()
        test.assert_equal(diff(cmn.Vec(0,0,0), cmn.Vec(1,0,0)), cmn.Vec(1,0,0),
                          err_msg="both vectors were already inside the domain")
        test.assert_equal(diff(cmn.Vec(0,0,0), cmn.Vec(-1.5, 0, 0)), cmn.Vec(.5, 0, 0),
                          err_msg="The latter vector was outside the domain, thus had to get a "
                                  "position update to (.5, 0, 0)")

    def test_kernel_context_dist_squared_fun(self):
        self.ctx.box_size = cmn.Vec(2,2,2)
        self.ctx.periodic_boundary = [True, True, True]
        dist = self.ctx.get_dist_squared_fun()
        test.assert_equal(dist(cmn.Vec(0,0,0), cmn.Vec(1,0,0)), 1.0)
        test.assert_equal(dist(cmn.Vec(-.5,0,0), cmn.Vec(-1.5, 0, 0)), 1.0)

    def test_kernel_context_particle_radius(self):
        self.ctx.set_diffusion_constant("A", 1.0)
        self.ctx.set_particle_radius("A", 5.0)
        test.assert_equal(self.ctx.get_particle_radius("A"), 5.0)

    def test_kernel_context_time_step(self):
        self.ctx.timestep = 111.0
        test.assert_equal(self.ctx.timestep, 111.0)

if __name__ == '__main__':
    unittest.main()
