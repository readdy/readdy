import unittest
import numpy as np

import readdy._internal.common as cmn
import readdy._internal.prototyping as pr
import readdy._internal.api as sim


class TestPrograms(unittest.TestCase):
    def test_custom_state_model(self):
        class CustomStateModel(pr.Model):
            def get_particle_positions(self):
                result = cmn.StdVectorVec3()
                result.append(cmn.Vec(-1, -1, -1))
                return result

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
        class CustomProgram(pr.Program):
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
        kernel.get_kernel_context().set_diffusion_constant("A", 1.0)
        factory = kernel.get_program_factory()

        add_particles = factory.create_add_particles()
        integrator = factory.create_euler_integrator()
        forces = factory.create_update_forces()
        neighbor_list = factory.create_update_neighbor_list()
        reactions = factory.create_reactions_uncontrolled_approximation()

        add_particles.add_particle(pr.Particle(0, 0, 0, 0))
        add_particles.add_particle(pr.Particle(1, 1, 1, 0))
        add_particles.execute()

        state_model = kernel.get_kernel_state_model()
        positions = state_model.get_particle_positions()
        np.testing.assert_equal(positions[0], cmn.Vec(0, 0, 0))
        np.testing.assert_equal(positions[1], cmn.Vec(1, 1, 1))

        it = kernel.get_kernel_state_model().get_particle_data().positions
        assert next(it) == cmn.Vec(0, 0, 0)
        assert next(it) == cmn.Vec(1, 1, 1)

        with np.testing.assert_raises(StopIteration):
            next(it)


if __name__ == '__main__':
    unittest.main()
