import unittest

import numpy as np
from readdy._internal.simulation import Vec
from readdy._internal.simulation import Simulation
from readdy._internal.simulation import KernelProvider

from readdy.util import platform_utils


class TestInternalSimulationModule(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.kernel_provider = KernelProvider.get()
        cls.kernel_provider.load_from_dir(platform_utils.get_readdy_plugin_dir())
        cls.simulation = Simulation()

    def test_properties(self):
        np.testing.assert_equal(self.simulation.isKernelSelected(), False)
        self.simulation.setKernel('SingleCPU')
        np.testing.assert_equal(self.simulation.isKernelSelected(), True)
        np.testing.assert_equal(self.simulation.getSelectedKernelType(), "SingleCPU")
        self.simulation.kbt = 5.0
        np.testing.assert_equal(self.simulation.kbt, 5.0)
        self.simulation.periodic_boundary = [True, False, True]
        np.testing.assert_equal(self.simulation.periodic_boundary, (True, False, True))
        self.simulation.box_size = Vec(1, 3.6, 7)
        np.testing.assert_equal(self.simulation.box_size, Vec(1, 3.6, 7))


if __name__ == '__main__':
    unittest.main()