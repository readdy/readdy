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

import readdy._internal.common as cmn
import readdy._internal.prototyping as pr
import readdy._internal.api as sim


class TestPrograms(unittest.TestCase):
    def test_custom_state_model(self):
        class CustomStateModel(pr.Model):
            def get_particle_positions(self):
                return [cmn.Vec(-1, -1, -1)]

        class CustomKernel(pr.SingleCPUKernel):
            def __init__(self):
                super(CustomKernel, self).__init__()
                self._model = CustomStateModel(self.get_kernel_context())

            def get_kernel_state_model(self):
                return self._model

        kernel = CustomKernel()
        pos = kernel.get_kernel_state_model().get_particle_positions()
        assert len(pos) == 1
        assert pos[0] == cmn.Vec(-1, -1, -1)

    def test_custom_program(self):
        class CustomProgram(pr.Action):
            def __init__(self, kernel):
                self._executed = False

            def execute(self):
                self._executed = True

            @property
            def executed(self):
                return self._executed

        kernel = pr.SingleCPUKernel()
        program = CustomProgram(kernel)
        program.execute()
        assert program.executed

    def test_factory_programs(self):
        kernel = pr.SingleCPUKernel()
        kernel.get_kernel_context().register_particle_type("A", 1.0, 1.0)
        factory = kernel.get_action_factory()

        add_particles = factory.create_add_particles([pr.Particle(0, 0, 0, 0), pr.Particle(1, 1, 1, 0)])
        integrator = factory.create_euler_integrator(1)
        forces = factory.create_update_forces()
        neighbor_list = factory.create_update_neighbor_list()
        reactions = factory.create_reactions_uncontrolled_approximation(1)

        add_particles.perform()

        state_model = kernel.get_kernel_state_model()
        positions = state_model.get_particle_positions()
        np.testing.assert_equal(positions[0], cmn.Vec(0, 0, 0))
        np.testing.assert_equal(positions[1], cmn.Vec(1, 1, 1))

        it = kernel.get_kernel_state_model().get_particle_data().entries
        assert next(it).pos == cmn.Vec(0, 0, 0)
        assert next(it).pos == cmn.Vec(1, 1, 1)

        with np.testing.assert_raises(StopIteration):
            next(it)


if __name__ == '__main__':
    unittest.main()
