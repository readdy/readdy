import os
import unittest

import numpy as np
import readdy2._internal.simulation


class TestInternalSimulationModule(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        maybe_prefix = os.environ.get('PREFIX')
        cls.kernel_provider = readdy2._internal.simulation.KernelProvider.get()
        if maybe_prefix is not None:
            # we are in the conda test environment
            plugin_dir = os.path.join(maybe_prefix, 'lib', 'readdy_plugins')
            cls.kernel_provider.load_from_dir(plugin_dir)
        cls.simulation = readdy2._internal.simulation.Simulation()
        cls.simulation.setKernel('SingleCPU')

    def test_properties(self):
        self.simulation.kbt = 5.0
        np.testing.assert_equal(self.simulation.kbt, 5.0)
        # todo test more properties
