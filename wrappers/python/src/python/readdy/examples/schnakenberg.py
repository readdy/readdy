# coding=utf-8

# Copyright © 2016 Computational Molecular Biology Group,
#                  Freie Universität Berlin (GER)
#
# This file is part of ReaDDy.
#
# ReaDDy is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General
# Public License along with this program. If not, see
# <http://www.gnu.org/licenses/>.


# coding=utf-8
from __future__ import print_function

"""
Erban, R., Chapman, J., & Maini, P. (2007). A practical guide to stochastic simulations of reaction-diffusion processes.
arXiv Preprint arXiv:0704.1908, 35. Retrieved from http://arxiv.org/abs/0704.1908

Page 30
"""

import numpy as np
from readdy._internal.api import KernelProvider, Simulation
from readdy._internal.common import Vec
from readdy.util import platform_utils

from scipy.optimize import brentq

import matplotlib.pyplot as plt


def erban_chapman(k, D, R, kappa):
    """
    Erban, Chapman:
    k - known rate in [x³/s]
    D := D_A + D_B (sum of diffusion coefficients)
    ĸ - microscopic rate [1/s]
    k = 4 * pi * D * (R - sqrt(D/ĸ)tanh(sqrt(ĸ/D)*R)
    """
    return k - 4 * np.pi * D * (R - np.sqrt(D / kappa) * np.tanh(np.sqrt(kappa / D) * R))


class SchnakenbergSimulation(object):
    def __init__(self, timesteps=100000):
        self._timesteps = timesteps
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111)
        self.ax.set_xscale('log')
        self.lines, = self.ax.plot([],[], '-')
        self.fig.show()
        plt.draw()
        plt.ioff()
        plt.pause(.01)

        self._data = np.zeros((timesteps+1, 2))
        self._t = 0

    def run(self):

        ###################################
        #
        # Units:
        #   - [t] = sec
        #
        ###################################

        tau = 1e-3

        kernel_provider = KernelProvider.get()
        kernel_provider.load_from_dir(platform_utils.get_readdy_plugin_dir())
        simulation = Simulation()
        simulation.set_kernel("CPU")
        simulation.kbt = 1.0
        simulation.box_size = Vec(3, 3, 3)
        simulation.periodic_boundary = [True, True, True]

        D = 1.0
        R = .01
        simulation.register_particle_type("A", D, R)
        simulation.register_particle_type("2A", D, R)
        simulation.register_particle_type("3A", D, R)
        simulation.register_particle_type("B", D, R)
        simulation.register_particle_type("GA", D, R)
        simulation.register_particle_type("GB", D, R)

        reaction_radius = .2
        k1 = 4 * 1e-5
        k2 = 50
        k3 = 10
        k4 = 250

        V = simulation.box_size * simulation.box_size

        N_GA = 100.0
        N_GB = 1000.0

        k_enzymatic = brentq(lambda x: erban_chapman(k1, 2, reaction_radius, x), 1e-10, 5000000000000)
        k_enzymatic = 10. # k_enzymatic for reaction_radius = .1
        print("k_enzymatic=%s" % k_enzymatic)
        print("2 * k3 - k3 * k3 * tau = %s" % (2.0 * k3 - k3 * k3 * tau))
        print("k3 * k3 * tau = %s" % (k3 * k3 * tau))
        print("k_birthA = %s" % (k2 * V / N_GA))
        print("k_birthB = %s" % (k4 * V / N_GB))
        print("sqrt(R*R/D) = %s" % (np.sqrt(reaction_radius * reaction_radius / D)))

        simulation.register_reaction_conversion("2A -> A", "2A", "A", 2.0 * k3 - k3 * k3 * tau)
        simulation.register_reaction_decay("2A -> 0", "2A", k3 * k3 * tau)
        simulation.register_reaction_fusion("A + A -> 2A", "A", "A", "2A", 1.0 / tau, reaction_radius, .5, .5)
        simulation.register_reaction_enzymatic("2A + B -> 2A + A", "2A", "B", "A", k_enzymatic, reaction_radius)
        simulation.register_reaction_fission("GA -> GA + A", "GA", "GA", "A", k2 * V / N_GA, reaction_radius, .5, .5)
        simulation.register_reaction_decay("A -> 0", "A", k3)
        simulation.register_reaction_fission("GB -> GB + B", "GB", "GB", "B", k4 * V / N_GB, reaction_radius, .5, .5)

        simulation.add_particle("A", Vec(0, 0, 0))

        ga_x = np.random.uniform(-0.5 * simulation.box_size[0], 0.5 * simulation.box_size[0], int(N_GA))
        ga_y = np.random.uniform(-0.5 * simulation.box_size[1], 0.5 * simulation.box_size[1], int(N_GA))
        ga_z = np.random.uniform(-0.5 * simulation.box_size[2], 0.5 * simulation.box_size[2], int(N_GA))

        gb_x = np.random.uniform(-0.5 * simulation.box_size[0], 0.5 * simulation.box_size[0], int(N_GB))
        gb_y = np.random.uniform(-0.5 * simulation.box_size[1], 0.5 * simulation.box_size[1], int(N_GB))
        gb_z = np.random.uniform(-0.5 * simulation.box_size[2], 0.5 * simulation.box_size[2], int(N_GB))

        for i in range(int(N_GA)):
            simulation.add_particle("GA", Vec(ga_x[i], ga_y[i], ga_z[i]))
        for i in range(int(N_GB)):
            simulation.add_particle("GB", Vec(gb_x[i], gb_y[i], gb_z[i]))

        N_B = 30000
        boxsize = np.array([simulation.box_size[0], simulation.box_size[1], simulation.box_size[2]])
        for i in range(N_B):
            pos = np.random.random(3) * boxsize - 0.5 * boxsize
            simulation.add_particle("B", Vec(pos[0], pos[1], pos[2]))

        def callback(result):
            n_a, n_b = result[0] + 2 * result[1] + 3 * result[2], result[3]
            self._data[self._t, 0] = n_a
            self._data[self._t, 1] = n_b
            print("(%s,%s,a=%s,a2=%s,a3=%s)"%(n_a, n_b,result[0], result[1], result[2]))

            self.lines.set_xdata(self._data[:self._t, 0])
            self.lines.set_ydata(self._data[:self._t, 1])
            self.ax.relim()
            self.ax.autoscale_view()
            #We need to draw *and* flush
            self.fig.canvas.draw()
            self.fig.canvas.flush_events()

            plt.pause(.1)
            self._t += 1

        simulation.register_observable_n_particles_types(1, ["A", "2A", "3A", "B"], callback)
        print("start run")
        simulation.run(self._timesteps, tau)


if __name__ == '__main__':
    simulation = SchnakenbergSimulation()
    simulation.run()
