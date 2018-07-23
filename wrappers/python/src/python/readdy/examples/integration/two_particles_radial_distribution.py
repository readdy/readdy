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


def plot_boltzmann(force_const, interaction_radius):
    def potential(r, force_const, interaction_radius):
        if r < interaction_radius:
            return 0.5 * force_const * np.power(r - interaction_radius, 2.)
        else:
            return 0.

    boltz = lambda r: np.exp(-1. * potential(r, force_const, interaction_radius))
    r_range = np.linspace(0.1, 2., 100)
    b_range = np.fromiter(map(boltz, r_range), dtype=float)
    plt.plot(r_range, b_range, label=r"Boltzmann correlation $e^{-\beta U(r)}$")


if __name__ == '__main__':
    system = readdy.ReactionDiffusionSystem(
        box_size=(4, 4, 4), periodic_boundary_conditions=[True, True, True], unit_system=None)
    system.add_species("A", 1.0)
    system.add_species("B", 1.0)
    system.potentials.add_harmonic_repulsion("A", "B", 1., 1.)

    simulation = system.simulation(kernel="SingleCPU")

    simulation.output_file = "out.h5"
    simulation.observe.rdf(200, np.linspace(0., 2., 10), ["A"], ["B"], 1. / system.box_volume)
    simulation.add_particle("A", [0., 0., 0.])
    simulation.add_particle("B", [0., 0., 1.])

    if os.path.exists(simulation.output_file):
        os.remove(simulation.output_file)

    simulation.run(n_steps=10000000, timestep=2e-3)

    traj = readdy.Trajectory(simulation.output_file)
    rdf_times, bin_centers, rdf_values = traj.read_observable_rdf()

    mean, std_dev, std_err = average_across_first_axis(rdf_values)
    plt.errorbar(bin_centers, mean, yerr=std_err, fmt=".", label="ReaDDy")
    plot_boltzmann(1., 1.)
    plt.legend()
    plt.xlabel(r"Distance $r$ of A and B")
    plt.ylabel(r"Radial distribution function $g(r)$")
    plt.show()
