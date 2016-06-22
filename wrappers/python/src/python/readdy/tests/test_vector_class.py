import unittest

from readdy._internal.simulation import Vec


class TestVectorClass(unittest.TestCase):

    def test_sanity(self):
        v = Vec(3, 3, 3)
        print (v * Vec(6, 6, 6))


if __name__ == '__main__':
    unittest.main()
