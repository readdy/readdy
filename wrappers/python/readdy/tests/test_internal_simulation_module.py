import os
import unittest

import numpy as np
import readdy._internal.simulation

from readdy.util import platform_utils


class TestInternalSimulationModule(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.kernel_provider = readdy._internal.simulation.KernelProvider.get()
        cls.kernel_provider.load_from_dir(platform_utils.get_readdy_plugin_dir())
        cls.simulation = readdy._internal.simulation.Simulation()
        cls.simulation.setKernel('SingleCPU')

    def test_properties(self):
        self.simulation.kbt = 5.0
        np.testing.assert_equal(self.simulation.kbt, 5.0)
        # todo test more properties

if __name__ == '__main__':
    unittest.main()