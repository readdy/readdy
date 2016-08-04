import unittest

import readdy._internal.prototyping as pr


class SingleCPUExtension(pr.SingleCPUKernel):
    def __init__(self):
        print("calling super init")
        super(SingleCPUExtension, self).__init__()

    def get_kernel_state_model(self):
        print("calling super get kernel state model")
        super(SingleCPUExtension, self).get_kernel_state_model()

def creator():
    return SingleCPUExtension()


class TestTest(unittest.TestCase):

    def test_load_prototyping_module(self):
        scpu = SingleCPUExtension()
        program_factory = scpu.get_program_factory()
        addp = program_factory.create_add_particles()
        addp.execute()


if __name__ == '__main__':
    unittest.main()