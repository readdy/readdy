# coding=utf-8

# Copyright © 2018 Computational Molecular Biology Group,
#                  Freie Universität Berlin (GER)
#
# Redistribution and use in source and binary forms, with or
# without modification, are permitted provided that the
# following conditions are met:
#  1. Redistributions of source code must retain the above
#     copyright notice, this list of conditions and the
#     following disclaimer.
#  2. Redistributions in binary form must reproduce the above
#     copyright notice, this list of conditions and the following
#     disclaimer in the documentation and/or other materials
#     provided with the distribution.
#  3. Neither the name of the copyright holder nor the names of
#     its contributors may be used to endorse or promote products
#     derived from this software without specific
#     prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
# CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
# INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
# STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
# ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

from __future__ import print_function
from readdy._internal.readdybinding.api import KernelProvider, Simulation
from readdy._internal.readdybinding.common import Vec

from readdy.util import platform_utils
from scipy.optimize import brentq

import numpy as np
from time import clock, time
import matplotlib

# matplotlib.use('Qt4Agg')

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


# abbreviation for plotting functions in certain range
def plot(f, x_range, **kw):
    y_range = np.array(map(f, x_range))
    plt.plot(x_range, y_range, **kw)


class MinEMinDSimulation(object):
    def __init__(self, result_fname=None, generate_plots=True):
        print("Init MinE/MinD simulation, fname=%s"%result_fname)
        self._result_fname = result_fname
        self._generate_plots = generate_plots
        self.t_d = 0
        self.t_e = 0
        self.stride = 2000
        self.timestep = .0005
        self.n_timesteps = 3000000
        self.time = None
        if self._generate_plots:
            self.fig = plt.figure()
            if(self._result_fname is not None):
                self.fig.suptitle(self._result_fname)
            self.axis = [self.fig.add_subplot(231), self.fig.add_subplot(232), self.fig.add_subplot(233), self.fig.add_subplot(234), self.fig.add_subplot(235), self.fig.add_subplot(236)]
            self.axis[0].set_title("MinD")
            self.axis[1].set_title("MinDP")
            self.axis[2].set_title("MinDPB")
            self.axis[3].set_title("MinE")
            self.axis[4].set_title("MinDE")
            self.axis[5].set_title("MinD+MinDP+MinDPB+MinDE")
            self.prev_pos = {}
            self.prev_mine = None
            self.prev_mind = None
            self.prev_minde = None
            self.fig.show()
            plt.draw()
            plt.ioff()
            plt.pause(.01)
        self._hist_data = [None, None, None, None, None, None]

    def callback_histogram(self, data, idx):

        if idx == 0:
            if self.time is not None:
                print("elapsed: {0}".format(time() - self.time))
            self.time = time()
        if self._hist_data[idx] is None:
            self._hist_data[idx] = np.array(data)
        else:
            self._hist_data[idx] = np.vstack((self._hist_data[idx], data))
            if self._generate_plots:
                self.axis[idx].imshow(self._hist_data[idx], cmap='hot')
                plt.pause(.1)
        if idx == 0:
            print("t={0} ({1} sec) -> {2:.3f}%".format(self.t_d*self.stride, self.t_d*self.stride*self.timestep, 100.*self.t_d*self.stride/float(self.n_timesteps)))
            self.t_d += 1

    def histogram_callback_minD(self, histogramTuple):
        self.callback_histogram(histogramTuple, 0)

    def histogram_callback_minDP(self, histogramTuple):
        self.callback_histogram(histogramTuple, 1)

    def histogram_callback_minDPB(self, histogramTuple):
        self.callback_histogram(histogramTuple, 2)

    def histogram_callback_minE(self, histogramTuple):
        self.callback_histogram(histogramTuple, 3)

    def histogram_callback_minDE(self, histogramTuple):
        self.callback_histogram(histogramTuple, 4)

    def histogram_callback_M(self, histogramTuple):
        self.callback_histogram(histogramTuple, 5)

    def n_particles_callback(self, n_particles):
        # ["D", "D_P", "D_PB", "E", "DE"]
        print("n_minD={0}, n_minDP={1}, n_minDPB={2}, n_minE={3}, n_minDE={4}, total={5}".format(n_particles[0], n_particles[1], n_particles[2], n_particles[3], n_particles[4], sum(n_particles)))

    def histrogram_callback_bound(self, histogramTuple):
        counts = histogramTuple[:]
        if self._hist_data[1] is None:
            self._hist_data[1] = np.array(counts)
        else:
            self._hist_data[1] = np.vstack((self._hist_data[1], counts))
            self.axis[2].imshow(self._hist_data[1], cmap='hot')
            plt.pause(.00001)

    # def com_callback_mind(self, centerOfMass):
    #     self.t_d += 1
    #     # print("center of mass = %s" % centerOfMass)
    #     if self.prev_mind is not None:
    #         xs, ys, zs = [self.prev_mind[0], centerOfMass[0]], [self.prev_mind[1], centerOfMass[1]], [self.prev_mind[2], centerOfMass[2]]
    #         self.ax.plot(xs, ys, zs, color="red")
    #         plt.pause(.001)
    #     self.prev_mind = centerOfMass
    #
    # def com_callback_mine(self, centerOfMass):
    #     self.t_e += 1
    #     # print("center of mass minE = %s" % centerOfMass)
    #     if self.prev_mine is not None:
    #         xs, ys, zs = [self.prev_mine[0], centerOfMass[0]], [self.prev_mine[1], centerOfMass[1]], [self.prev_mine[2], centerOfMass[2]]
    #         self.ax.plot(xs, ys, zs, color="blue")
    #         plt.pause(.001)
    #     self.prev_mine = centerOfMass
    #
    # def com_callback_minde(self, centerOfMass):
    #     if self.prev_minde is not None:
    #         xs, ys, zs = [self.prev_minde[0], centerOfMass[0]], [self.prev_minde[1], centerOfMass[1]], [self.prev_minde[2], centerOfMass[2]]
    #         self.ax.plot(xs, ys, zs, color="orange")
    #         plt.pause(.001)
    #     self.prev_minde = centerOfMass

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
        simulation.set_kernel("CPU")

        ###################################
        #
        # set up simulation box
        #
        ###################################

        box_size = Vec(2, 7, 12)
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

        membrane_particle_size = .05

        diffusion_factor = .5
        simulation.register_particle_type("D", 2.5 * diffusion_factor, .01)  # MinD-ADP (without phosphor)
        simulation.register_particle_type("D_P", 2.5 * diffusion_factor, .01)  # MinD-ATP (with phosphor)
        simulation.register_particle_type("E", 2.5 * diffusion_factor, .01)  # MinE
        simulation.register_particle_type("D_PB", .01 * diffusion_factor, .01)  # MinD-ATP bound
        simulation.register_particle_type("DE", .01 * diffusion_factor, .01)  # MinDE

        ###################################
        #
        # register reaction types
        #
        ###################################

        reaction_radius = 4*(0.01 + 0.01)  # = sum of the particle radii * 5 (5 - magic number such that k_fusion makes sense, sort of) 5 *
        # k_fusion = brentq(lambda x: self.erban_chapman(.093, 2.5 + .01, reaction_radius, x), 1, 5000000)
        k_fusion = 1.0
        print("k_fusion=%s" % k_fusion)
        simulation.register_reaction_conversion("Phosphorylation", "D", "D_P", .5)
        simulation.register_reaction_fusion("bound MinD+MinE->MinDE", "D_PB", "E", "DE", k_fusion, reaction_radius*3.5, .5, .5)
        simulation.register_reaction_fission("MinDE to MinD and MinE, detach", "DE", "D", "E", .25, reaction_radius, .5, .5)

        ###################################
        #
        # register potentials
        #
        ###################################

        membrane_size = Vec(.5, 5, 10)
        layer = Vec(.08, .08, .08)
        extent = membrane_size + 2 * layer
        origin = -.5 * membrane_size - layer
        simulation.register_potential_box("D", 10., origin, extent, False)  # (force constant, origin, extent, considerParticleRadius)
        simulation.register_potential_box("D_P", 10., origin, extent, False)  # (force constant, origin, extent, considerParticleRadius)
        simulation.register_potential_box("D_PB", 10., origin, extent, False)  # (force constant, origin, extent, considerParticleRadius)
        simulation.register_potential_box("E", 10., origin, extent, False)  # (force constant, origin, extent, considerParticleRadius)
        simulation.register_potential_box("DE", 10., origin, extent, False)  # (force constant, origin, extent, considerParticleRadius)

        # simulation.register_potential_piecewise_weak_interaction("D_P", "D_PB", 3, .02, 2, .05)  # (force constant, desired dist, depth, no interaction dist)

        ###################################
        #
        # membrane particles
        #
        ###################################
        using_membrane_particles = False
        if using_membrane_particles:
            simulation.register_particle_type("M", 0, membrane_particle_size)  # membrane particle
            simulation.register_reaction_enzymatic("Attach to membrane", "M", "D_P", "D_PB", .5, .01 + membrane_particle_size)  # .01 + .025  # todo: rate?
            dx = np.linspace(origin[0] + layer[0], -1 * origin[0] - layer[0], int(float(membrane_size[0]) / membrane_particle_size), endpoint=True)
            dy = np.linspace(origin[1] + layer[1], -1 * origin[1] - layer[1], int(float(membrane_size[1]) / membrane_particle_size), endpoint=True)
            dz = np.linspace(origin[2] + layer[2], -1 * origin[2] - layer[2], int(float(membrane_size[2]) / membrane_particle_size), endpoint=True)
            for y in dy:
                for z in dz:
                    simulation.add_particle("M", Vec(-1 * origin[0] - layer[0], y, z))
            print("done adding membrane particles")
        else:
            simulation.register_reaction_conversion("Phosphorylation", "D_P", "D_PB", .5)
            simulation.register_reaction_enzymatic("Enzymatic DP+DPB->DPB + DPB", "D_PB", "D_P", "D_PB", .5, .02)
        using_uniform_distribution = True
        n_minE_particles = 3120
        n_minD_particles = n_minE_particles * 4
        mine_x = np.random.uniform(origin[0] + layer[0], -1 * origin[0] - layer[0], n_minE_particles)
        mine_y = np.random.uniform(origin[1] + layer[1], -1 * origin[1] - layer[1], n_minE_particles)
        if using_uniform_distribution:
            mine_z = np.random.uniform(origin[2] + layer[2], -1 * origin[2] - layer[2], n_minE_particles)
        else:
            mine_z = np.random.uniform(origin[2] + layer[2], .5 * (-1 * origin[2] - layer[2]), n_minE_particles)

        mind_x = np.random.uniform(origin[0] + layer[0], -1 * origin[0] - layer[0], n_minD_particles)
        mind_y = np.random.uniform(origin[1] + layer[1], -1 * origin[1] - layer[1], n_minD_particles)
        if using_uniform_distribution:
            mind_z = np.random.uniform(origin[2] + layer[2], -1 * origin[2] - layer[2], n_minD_particles)
        else:
            mind_z = np.random.uniform(.5 * (-1 * origin[2] - layer[2]), -1 * origin[2] - layer[2], n_minD_particles)

        for i in range(n_minE_particles):
            simulation.add_particle("E", Vec(mine_x[i], mine_y[i], mine_z[i]))

        for i in range(int(.5 * n_minD_particles)):
            simulation.add_particle("D", Vec(mind_x[i], mind_y[i], mind_z[i]))
        for i in range(int(.5 * n_minD_particles), n_minD_particles):
            simulation.add_particle("D_P", Vec(mind_x[i], mind_y[i], mind_z[i]))

        self.timestep = simulation.get_recommended_time_step(2)

        ###################################
        #
        # register observables
        #
        ###################################

        # simulation.register_observable_center_of_mass(1, self.com_callback_mind, ["D", "D_P", "D_PB"])
        # simulation.register_observable_center_of_mass(1, self.com_callback_mine, ["E"])
        # simulation.register_observable_center_of_mass(1, self.com_callback_minde, ["DE", "D_PB"])
        print("histogram start")
        # simulation.register_observable_histogram_along_axis(100, self.histrogram_callback_minD, np.arange(-3, 3, .1), ["D", "D_P", "D_PB"], 2)
        # simulation.register_observable_histogram_along_axis(100, self.histrogram_callback_minE, np.arange(-3, 3, .1), ["D_PB", "DE"], 2)
        stride = int(.01/self.timestep)
        self.stride = stride
        print("using stride=%s" % stride)
        bins = np.linspace(-7, 7, 80)
        simulation.register_observable_histogram_along_axis(stride, bins, 2, ["D"], self.histogram_callback_minD)
        simulation.register_observable_histogram_along_axis(stride, bins, 2, ["D_P"], self.histogram_callback_minDP)
        simulation.register_observable_histogram_along_axis(stride, bins, 2, ["D_PB"], self.histogram_callback_minDPB)
        simulation.register_observable_histogram_along_axis(stride, bins, 2, ["E"], self.histogram_callback_minE)
        simulation.register_observable_histogram_along_axis(stride, bins, 2, ["DE"], self.histogram_callback_minDE)
        simulation.register_observable_histogram_along_axis(stride, bins, 2, ["D", "D_P", "D_PB", "DE"], self.histogram_callback_M)
        simulation.register_observable_n_particles(stride, ["D", "D_P", "D_PB", "E", "DE"], self.n_particles_callback)
        print("histogram end")

        self.n_timesteps = int(1200./self.timestep)

        print("starting simulation for effectively %s sec" % (self.timestep * self.n_timesteps))
        simulation.run_scheme_readdy(True).with_reaction_scheduler("GillespieParallel").configure(self.timestep).run(self.n_timesteps)

        if self._result_fname is not None:
            with open(self._result_fname, 'w') as f:
                np.save(f, np.array(self._hist_data))


if __name__ == '__main__':
    sim = MinEMinDSimulation('test_mind_mine.npy', True)
    sim.execute()

    #X = np.load('test_mind_mine_no_membrane6.npy')
    #fig = plt.figure()
    #ax = fig.add_subplot(111)
    #print(X[1].shape)
    #ax.imshow(X[5][450:], cmap='hot')
    #plt.show()
