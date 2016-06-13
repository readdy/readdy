import unittest

from readdy._internal.simulation import Vec


class TestVectorClass(unittest.TestCase):
    def calcEnergy(self, v1, v2):
        print("calc energy with v1=%s, v2=%s" % (v1, v2))
        return 0

    def calcForce(self, v1, v2):
        return Vec(0, 0, 0)

    def test_sanity(self):
        v = Vec(3, 3, 3)
        print (v * Vec(6, 6, 6))


if __name__ == '__main__':
    unittest.main()
