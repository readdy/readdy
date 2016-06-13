import unittest

from readdy._internal.simulation import Vec
from readdy._internal.simulation import Pot2


class TestVectorClass(unittest.TestCase):
    def calcEnergy(self, v1, v2):
        print("calc energy with v1=%s, v2=%s" % (v1, v2))
        return 0

    def calcForce(self, v1, v2):
        return Vec(0, 0, 0)

    def test_test(self):
        pot = Pot2("mypot", self.calcEnergy, self.calcForce)

        pot.calc_energy(Vec(5, 5, 5), Vec(3, 3, 3))

    def test_sanity(self):
        v = Vec(3, 3, 3)
        print (v * Vec(6, 6, 6))


if __name__ == '__main__':
    unittest.main()
