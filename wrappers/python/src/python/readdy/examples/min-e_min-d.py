# coding=utf-8
from __future__ import print_function
from readdy._internal.simulation import KernelProvider, Simulation, Vec

from readdy.util import platform_utils
from scipy.optimize import brentq

import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


# abbreviation for plotting functions in certain range
def plot(f, x_range, **kw):
    y_range = np.array(map(f, x_range))
    plt.plot(x_range, y_range, **kw)


class MinEMinDSimulation(object):
    def __init__(self):
        self.t_d = 0
        self.t_e = 0
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111, projection='3d')
        plt.xlabel("x")
        plt.ylabel("y")
        self.prev_pos = {}
        self.fig.show()
        self.prev_mine = None
        self.prev_mind = None
        self.prev_minde = None
        plt.ioff()

    def com_callback_mind(self, centerOfMass):
        self.t_d += 1
        # print("center of mass = %s" % centerOfMass)
        if self.prev_mind is not None:
            xs, ys, zs = [self.prev_mind[0], centerOfMass[0]], [self.prev_mind[1], centerOfMass[1]], [self.prev_mind[2], centerOfMass[2]]
            self.ax.plot(xs, ys, zs, color="red")
            plt.pause(.001)
        self.prev_mind = centerOfMass

    def com_callback_mine(self, centerOfMass):
        self.t_e += 1
        # print("center of mass minE = %s" % centerOfMass)
        if self.prev_mine is not None:
            xs, ys, zs = [self.prev_mine[0], centerOfMass[0]], [self.prev_mine[1], centerOfMass[1]], [self.prev_mine[2], centerOfMass[2]]
            self.ax.plot(xs, ys, zs, color="blue")
            plt.pause(.001)
        self.prev_mine = centerOfMass

    def com_callback_minde(self, centerOfMass):
        if self.prev_minde is not None:
            xs, ys, zs = [self.prev_minde[0], centerOfMass[0]], [self.prev_minde[1], centerOfMass[1]], [self.prev_minde[2], centerOfMass[2]]
            self.ax.plot(xs, ys, zs, color="orange")
            plt.pause(.001)
        self.prev_minde = centerOfMass

    def erban_chapman(self, k, D, R, kappa):
        """
        Erban, Chapman:
        k - known rate in [x³/s]
        D := D_A + D_B (sum of diffusion coefficients)
        ĸ - microscopic rate [1/s]
        k = 4 * pi * D * (R - sqrt(D/ĸ)tanh(sqrt(ĸ/D)*R)
        """
        return k - 4 * np.pi * D * (R - np.sqrt(D / kappa) * np.tanh(np.sqrt(kappa / D) * R))

    def execute(self):
        ###################################
        #
        # Units:
        #   - [x] = µm
        #   - [t] = s
        #   - [E] = kJ/mol
        #
        ###################################


        kernel_provider = KernelProvider.get()
        kernel_provider.load_from_dir(platform_utils.get_readdy_plugin_dir())
        simulation = Simulation()
        simulation.setKernel("SingleCPU")

        ###################################
        #
        # set up simulation box
        #
        ###################################

        box_size = Vec(10, 10, 20)
        simulation.box_size = box_size
        simulation.kbt = 2.437  # room temperature
        simulation.periodic_boundary = [False, False, False]

        ###################################
        #
        # register particle types
        #
        ###################################

        # particle size, see: http://bmccellbiol.biomedcentral.com/articles/10.1186/1471-2121-5-29
        # "The size of the V-ATPase complex is about 15 nm (diameter) x 25 nm (length from lumen side to tip of head)"

        simulation.registerParticleType("M", 0, .125)  # membrane particle
        simulation.registerParticleType("D", 2.5, .01)  # MinD-ADP (without phosphor)
        simulation.registerParticleType("D_P", 2.5, .01)  # MinD-ATP (with phosphor)
        simulation.registerParticleType("E", 2.5, .01)  # MinE
        simulation.registerParticleType("D_PB", .01, .01)  # MinD-ATP bound
        simulation.registerParticleType("DE", .01, .01)  # MinDE

        ###################################
        #
        # register reaction types
        #
        ###################################

        reaction_radius = 5 * (0.01 + 0.01)  # = sum of the particle radii * 5 (5 - magic number such that k_fusion makes sense, sort of)
        k_fusion = brentq(lambda x: self.erban_chapman(.093, 2.5 + .01, reaction_radius, x), 1, 50)
        print("k_fusion=%s" % k_fusion)
        simulation.registerConversionReaction("Phosphorylation", "D", "D_P", .5)
        simulation.registerEnzymaticReaction("Attach to membrane", "M", "D_P", "D_PB", .5, .01 + .025)  # todo: rate?
        simulation.registerFusionReaction("bound MinD+MinE->MinDE", "D_PB", "E", "DE", k_fusion, reaction_radius, .5, .5)
        simulation.registerFissionReaction("MinDE to MinD and MinE, detach", "DE", "D", "E", .25, reaction_radius, .5, .5)

        ###################################
        #
        # register observables
        #
        ###################################

        simulation.registerObservable_CenterOfMass(1, self.com_callback_mind, ["D", "D_P", "D_PB"])
        simulation.registerObservable_CenterOfMass(1, self.com_callback_mine, ["E"])
        simulation.registerObservable_CenterOfMass(1, self.com_callback_minde, ["DE"])

        ###################################
        #
        # register potentials
        #
        ###################################

        membrane_size = Vec(2, 2, 4)
        extent = membrane_size + Vec(1, 1, 1)
        origin = -.5 * membrane_size - Vec(.5, .5, .5)
        simulation.registerBoxPotential("D", 10., origin, extent, True)  # (force constant, origin, extent, considerParticleRadius)
        simulation.registerBoxPotential("D_P", 10., origin, extent, True)  # (force constant, origin, extent, considerParticleRadius)
        simulation.registerBoxPotential("D_PB", 10., origin, extent, True)  # (force constant, origin, extent, considerParticleRadius)
        simulation.registerBoxPotential("E", 10., origin, extent, True)  # (force constant, origin, extent, considerParticleRadius)
        simulation.registerBoxPotential("DE", 10., origin, extent, True)  # (force constant, origin, extent, considerParticleRadius)

        simulation.registerWeakInteractionPiecewiseHarmonicPotential("D_P", "D_PB", 10, .025, 2, .05)  # (force constant, desired dist, depth, no interaction dist)

        ###################################
        #
        # membrane particles
        #
        ###################################
        stepwidth = 0.2
        dx = dy = np.arange(origin[0] + .5, -1 * origin[0] - .5 + stepwidth, stepwidth)
        dz = np.arange(origin[2] + .5, -1 * origin[2] - .5 + stepwidth, stepwidth)
        # front and back
        for x in dx:
            for y in dy:
                simulation.addParticle("M", Vec(x, y, origin[2] + .5))
                simulation.addParticle("M", Vec(x, y, -1 * origin[2] - .5))
        for x in dx:
            for z in dz:
                simulation.addParticle("M", Vec(x, origin[1] + .5, z))
                simulation.addParticle("M", Vec(x, -1 * origin[1] - .5, z))
                simulation.addParticle("M", Vec(origin[0] + .5, x, z))
                simulation.addParticle("M", Vec(-1 * origin[0] - .5, x, z))
        print("done adding membrane particles")
        n_minE_particles = 50
        n_minD_particles = 50*4
        mine_x = np.random.uniform(origin[0] + .5, -1 * origin[0] - .5, n_minE_particles)
        mine_y = np.random.uniform(origin[1] + .5, -1 * origin[1] - .5, n_minE_particles)
        mine_z = np.random.uniform(origin[2] + .5, -1 * origin[2] - .5, n_minE_particles)

        mind_x = np.random.uniform(origin[0] + .5, -1 * origin[0] - .5, n_minD_particles)
        mind_y = np.random.uniform(origin[1] + .5, -1 * origin[1] - .5, n_minD_particles)
        mind_z = np.random.uniform(origin[2] + .5, -1 * origin[2] - .5, n_minD_particles)

        for i in range(n_minE_particles):
            simulation.addParticle("E", Vec(mine_x[i], mine_y[i], mine_z[i]))

        for i in range(n_minD_particles):
            simulation.addParticle("D", Vec(mind_x[i], mind_y[i], mind_z[i]))

        # membrane_positions = simulation.getParticlePositions("M")
        # print("done retrieving positions")
        # xs,ys,zs = np.array([v[0] for v in membrane_positions]), np.array([v[1] for v in membrane_positions]), np.array([v[2] for v in membrane_positions])
        # print("done transposing")
        # fig = plt.figure()
        # ax = fig.add_subplot(111, projection='3d')
        # ax.scatter(xs,ys,zs)
        # plt.show()

        print("starting simulation")
        simulation.run(10000, .05)


if __name__ == '__main__':
    sim = MinEMinDSimulation()
    sim.execute()
