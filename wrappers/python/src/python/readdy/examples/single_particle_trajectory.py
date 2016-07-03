import unittest

from readdy._internal.simulation import KernelProvider, Simulation, Vec

from readdy.util import platform_utils
import matplotlib

matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class TestPlot(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.kernel_provider = KernelProvider.get()
        cls.kernel_provider.load_from_dir(platform_utils.get_readdy_plugin_dir())
        cls.simulation = Simulation()
        cls.simulation.setKernel("SingleCPU")

        cls.fig = plt.figure()
        cls.ax = cls.fig.add_subplot(111, projection='3d')
        cls.prev_pos = {}
        cls.fig.show()
        cls.ax.set_xlim([-3, 3])
        cls.ax.set_ylim([-3, 3])
        cls.ax.set_zlim([-3, 3])
        cls.current_plot = None
        plt.ioff()

    def position_callback(self, positions):

        for key in self.prev_pos.keys():
            pp = self.prev_pos[key]
            xs, ys, zs = [pp[0], positions[key][0]], [pp[1], positions[key][1]], [pp[2], positions[key][2]]
            self.ax.plot(xs, ys, zs, color="blue")
        for i in range(len(positions)):
            self.prev_pos[i] = positions[i]

        plt.pause(.1)

    def test_plot(self):
        self.simulation.kbt = 1.0
        self.simulation.periodic_boundary = [False, False, False]
        self.simulation.box_size = Vec(5.0, 5.0, 5.0)
        self.simulation.registerParticleType("A", .1, 1.0)
        self.simulation.registerHarmonicRepulsionPotential("A", "A", .1)
        self.simulation.registerBoxPotential("A", 10, Vec(-1, -1, -1), Vec(2, 2, 2), True)
        self.simulation.addParticle("A", Vec(0, 0, 0))
        self.simulation.registerObservable_ParticlePositions(1, self.position_callback)
        self.simulation.run(100, .1)

if __name__ == '__main__':
    unittest.main()
