import unittest

from readdy._internal.api import Simulation

class TestSchemeApi(unittest.TestCase):

    def test_sanity(self):
        simulation = Simulation()
        simulation.set_kernel("SingleCPU")
        simulation.runReaDDyScheme(False)

if __name__ == '__main__':
    unittest.main()