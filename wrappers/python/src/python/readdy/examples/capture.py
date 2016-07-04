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

        self._rdf_centers = None
        self._rdf = None
        self._rdf_n_calls = 0

        self.T = 4000000

    def harmonic_force(self, r, forceConstant=5., depth=2., desired_dist=.2, no_interaction_dist=.5):
        len_part2 = no_interaction_dist - desired_dist
        if r < desired_dist:
            return -1 * forceConstant * (desired_dist - r)
        else:
            if r < desired_dist + .5 * len_part2:
                return -1 * depth * (1 / (.5 * len_part2)) * (1 / (.5 * len_part2)) * (r - desired_dist)
            else:
                if r < no_interaction_dist:
                    return depth * (1 / (.5 * len_part2)) * (1 / (.5 * len_part2)) * (no_interaction_dist - r)
        return 0

    def harmonic_energy(self, r, forceConstant=5., depth=2., desired_dist=.2, no_interaction_dist=.5):
        len_part2 = no_interaction_dist - desired_dist
        if r < desired_dist:
            return .5 * forceConstant * (r - desired_dist) * (r - desired_dist) - depth
        else:
            if r < desired_dist + .5 * len_part2:
                return .5 * depth * (1 / (.5 * len_part2)) * (1 / (.5 * len_part2)) * (r - desired_dist) * (r - desired_dist) - depth
            else:
                if r < no_interaction_dist:
                    return -.5 * depth * (1 / (.5 * len_part2)) * (1 / (.5 * len_part2)) * (r - no_interaction_dist) * (r - no_interaction_dist)
        return 0

    def ppos_callback(self, pos):
        if self._plot_particles:
            plt.cla()
            self.ax.set_xlim([-3, 3])
            self.ax.set_ylim([-3, 3])
            self.ax.set_zlim([-3, 3])
            r = [-.5, .5]
            for s, e in combinations(np.array(list(product(r, r, r))), 2):
                if np.sum(np.abs(s - e)) == r[1] - r[0]:
                    self.ax.plot3D(*zip(s, e), color="b")
            A = self.simulation.getParticlePositions("A")[0]
            B = self.simulation.getParticlePositions("B")[0]
            self.ax.scatter([A[0]], [A[1]], [A[2]], color="g", s=100)
            self.ax.scatter([B[0]], [B[1]], [B[2]], color="r", s=100)
            plt.pause(.001)

    def rdf_callback(self, pair):
        if self._rdf_centers is None:
            self._rdf_centers = pair[0]
        if self._rdf is None:
            self._rdf = pair[1]
        else:
            self._rdf += pair[1]
        self._rdf_n_calls += 1
        if self._rdf_n_calls % 10000 == 0:
            print("%s" % (10. * float(self._rdf_n_calls) / float(self.T)))

    def start(self):
        box_size = Vec(1.0, 1.0, 1.0)
        depth = .5
        desired_dist = .25
        force_constant = 4 * depth / (desired_dist * desired_dist)
        no_interaction_dist = 0.4
        print("fc=%s" % force_constant)
        self.simulation.kbt = 1.
        self.simulation.periodic_boundary = [True, True, True]
        self.simulation.box_size = box_size
        self.simulation.registerParticleType("A", .1, .1)
        self.simulation.registerParticleType("B", .1, .1)
        self.simulation.registerWeakInteractionPiecewiseHarmonicPotential("A", "B", force_constant, desired_dist, depth,
                                                                          no_interaction_dist)  # (force constant, desired dist, depth, no interaction dist)
        self.simulation.registerObservable_RadialDistribution(10, self.rdf_callback, np.arange(0, .5, .01), "A", "B", 1. / (box_size[0] * box_size[1] * box_size[2]))
        # self.simulation.registerBoxPotential("A", 10, Vec(-2, -2, -2), Vec(4, 4, 4), True)
        # self.simulation.registerBoxPotential("B", 10, Vec(-2, -2, -2), Vec(4, 4, 4), True)
        self.simulation.addParticle("A", Vec(-.0, -.0, -.0))
        self.simulation.addParticle("B", Vec(.5, .5, .5))
        self.simulation.registerObservable_ParticlePositions(1, self.ppos_callback)

        self.simulation.run(self.T, .001)

        x = np.arange(0, .5, .001)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(self._rdf_centers, self._rdf / self._rdf_n_calls)
        # ax.plot(x, map(lambda r: self.harmonic_energy(r, forceConstant=force_constant, depth=depth, desired_dist=desired_dist, no_interaction_dist=.3), x))
        ax.plot(x, map(lambda r: np.exp(-1 * self.harmonic_energy(r, forceConstant=force_constant, depth=depth, desired_dist=desired_dist, no_interaction_dist=no_interaction_dist)), x))
        # ax.plot(x, map(lambda r: self.harmonic_force(r, forceConstant=force_constant, depth=depth, desired_dist=desired_dist, no_interaction_dist=no_interaction_dist), x))
        plt.show()


if __name__ == '__main__':
    c = Capture(plot_particles=True)
    c.start()
