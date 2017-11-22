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


import os
import numpy as np
import matplotlib.pyplot as plt

import readdy


def average_across_first_axis(values):
    n_values = len(values)
    mean = np.sum(values, axis=0) / n_values  # shape = n_bins
    difference = values - mean  # broadcasting starts with last axis
    std_dev = np.sqrt(np.sum(difference * difference, axis=0) / n_values)
    std_err = np.sqrt(np.sum(difference * difference, axis=0) / n_values ** 2)
    return mean, std_dev, std_err


if __name__ == '__main__':
    system = readdy.ReactionDiffusionSystem(
        box_size=(10, 10, 10), periodic_boundary_conditions=[False, False, False], unit_system=None)
    system.add_species("A", 1.0)
    system.potentials.add_box(
        particle_type="A", force_constant=50., origin=[-4, -4, -4], extent=[8, 8, 8]
    )

    simulation = system.simulation(kernel="SingleCPU")

    simulation.output_file = "out.h5"
    simulation.observe.particle_positions(stride=1)
    init_pos = np.zeros((400, 3))
    simulation.add_particles("A", init_pos)

    if os.path.exists(simulation.output_file):
        os.remove(simulation.output_file)

    simulation.run(n_steps=5000, timestep=1e-3)

    traj = readdy.Trajectory(simulation.output_file)
    times, positions = traj.read_observable_particle_positions()
    times = np.array(times) * 1e-3

    # convert to pure numpy array to make use of fancy operations
    T = len(positions)
    N = len(positions[0])
    pos = np.zeros(shape=(T, N, 3))
    for t in range(T):
        for n in range(N):
            pos[t, n, 0] = positions[t][n][0]
            pos[t, n, 1] = positions[t][n][1]
            pos[t, n, 2] = positions[t][n][2]

    difference = pos - init_pos
    squared_displacements = np.sum(difference * difference, axis=2)  # sum over coordinates, per particle per timestep
    squared_displacements = squared_displacements.transpose()  # T x N -> N x T

    mean, std_dev, std_err = average_across_first_axis(squared_displacements)
    stride = 50
    plt.errorbar(times[::stride], mean[::stride], yerr=std_err[::stride], fmt=".", label="Particle diffusion")
    plt.plot(times[::stride], 6. * times[::stride], label=r"$6 D t$")
    plt.plot([stride * 1e-3, times[-1]], [16., 16.], "--", label="closest system boundary")
    plt.legend(loc="best")
    plt.xlabel(r"Time $t$")
    plt.ylabel(r"Mean squared displacement $\langle (x(t) - x_0)^2 \rangle_N$")
    plt.show()
