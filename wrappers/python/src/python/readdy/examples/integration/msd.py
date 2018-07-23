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
