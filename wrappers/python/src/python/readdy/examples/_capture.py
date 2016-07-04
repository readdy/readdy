from itertools import combinations, product

from readdy._internal.simulation import KernelProvider, Simulation, Vec

from readdy.util import platform_utils
import numpy as np
import matplotlib

matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class Capture(object):
    def __init__(self, plot_particles=False):
        self.kernel_provider = KernelProvider.get()
        self.kernel_provider.load_from_dir(platform_utils.get_readdy_plugin_dir())
        self.simulation = Simulation()
        self.simulation.setKernel("SingleCPU")
        self._plot_particles = plot_particles

        if self._plot_particles:
            self.fig = plt.figure()
            self.ax = self.fig.add_subplot(111, projection='3d')
            plt.ioff()
            self.fig.show()
            self.prev_pos = {}
            self.current_plot = None

        self.T = 4000000

    def ppos_callback(self, pos):
        if self._plot_particles:
            plt.cla()
            self.ax.set_xlim([-1, 1])
            self.ax.set_ylim([-1, 1])
            self.ax.set_zlim([-1, 1])
            r = [-.75, .75]
            for s, e in combinations(np.array(list(product(r, r, r))), 2):
                if np.sum(np.abs(s - e)) == r[1] - r[0]:
                    self.ax.plot3D(*zip(s, e), color="b")
            A = self.simulation.getParticlePositions("A")[0]
            B = self.simulation.getParticlePositions("B")[0]
            self.ax.scatter([A[0]], [A[1]], [A[2]], color="g", s=100)
            self.ax.scatter([B[0]], [B[1]], [B[2]], color="r", s=100)
            self.ax.plot3D([A[0], B[0]], [A[1], B[1]], [A[2], B[2]], color="r")
            plt.pause(.001)

    def start(self):
        box_size = Vec(2.0, 2.0, 2.0)
        depth = 2.
        desired_dist = .25
        force_constant = 4 * depth / (desired_dist * desired_dist)
        no_interaction_dist = 1.5
        print("fc=%s" % force_constant)
        self.simulation.kbt = 0.01
        self.simulation.periodic_boundary = [False, False, False]
        self.simulation.box_size = box_size
        self.simulation.registerParticleType("A", .1, .1)
        self.simulation.registerParticleType("B", .01, .1)
        self.simulation.registerWeakInteractionPiecewiseHarmonicPotential("A", "B", force_constant, desired_dist, depth,
                                                                          no_interaction_dist)  # (force constant, desired dist, depth, no interaction dist)
        self.simulation.registerBoxPotential("A", 100., Vec(-.75, -.75, -.75), Vec(1.5, 1.5, 1.5), False)
        self.simulation.registerBoxPotential("B", 100., Vec(-.75, -.75, -.75), Vec(1.5, 1.5, 1.5), False)
        self.simulation.addParticle("A", Vec(-.0, -.0, -.0))
        self.simulation.addParticle("B", Vec(0.1, 0.1, 0.1))
        self.simulation.registerObservable_ParticlePositions(1, self.ppos_callback)

        self.simulation.run(self.T, .0001)


if __name__ == '__main__':
    c = Capture(plot_particles=True)
    c.start()
