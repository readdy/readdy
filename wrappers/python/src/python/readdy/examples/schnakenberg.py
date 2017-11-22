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


"""
Erban, R., Chapman, J., & Maini, P. (2007). A practical guide to stochastic simulations of reaction-diffusion processes.
arXiv Preprint arXiv:0704.1908, 35. Retrieved from http://arxiv.org/abs/0704.1908

Page 30
"""

import os
import numpy as np
import matplotlib.pyplot as plt
# from scipy.optimize import brentq

import readdy


def k_erbchap(D, R, lambda_on):
    return 4 * np.pi * D * (R - np.sqrt(D / lambda_on) * np.tanh(np.sqrt(lambda_on / D) * R))


if __name__ == '__main__':
    system = readdy.ReactionDiffusionSystem(box_size=(10, 10, 10), periodic_boundary_conditions=[True, True, True],
                                            unit_system={"length_unit": "micrometer", "time_unit": "second"})
    system.add_species("A", 10.0)
    system.add_species("2A", 10.0)
    system.add_species("B", 10.0)
    system.add_species("FA", 10.0)
    system.add_species("FB", 10.0)
    n_factory_a = 50
    n_factory_b = 50
    lambda_on = 8.  # brentq(lambda x: 4e-2 - k_erbchap(3., 1., x), 1e-10, 5000000000000)
    production_rate_a = 50.
    decay_rate_a = 10.
    production_rate_b = 200.

    print("lambda_on", lambda_on)
    system.reactions.add("eat: B +(1) 2A -> A + 2A", rate=lambda_on)
    system.reactions.add("produce_a: FA -> FA +(1) A", rate=production_rate_a / n_factory_a)
    system.reactions.add("produce_b: FB -> FB +(1) B", rate=production_rate_b / n_factory_b)
    system.reactions.add("decay_a: A ->", rate=decay_rate_a)
    system.reactions.add("decay_2a: 2A -> A", rate=2. * decay_rate_a)
    system.reactions.add("decay_2a_2: 2A -> ", rate=1e-3 * decay_rate_a * decay_rate_a)
    system.reactions.add("fusion: A +(1) A -> 2A", rate=100000.)

    simulation = system.simulation(kernel="CPU")

    simulation.reaction_handler = "UncontrolledApproximation"
    simulation.output_file = "out.h5"
    simulation.observe.number_of_particles(stride=10, types=["A", "2A", "B", "FA", "FB"],
                                           callback=lambda c: print("B: ", c[2], " A: ", c[0] + 2 * c[1]))
    init_pos_fa = np.random.random(size=(n_factory_a, 3)) * 10. - 5.
    init_pos_fb = np.random.random(size=(n_factory_b, 3)) * 10. - 5.
    simulation.add_particles("FA", init_pos_fa)
    simulation.add_particles("FB", init_pos_fb)

    if os.path.exists(simulation.output_file):
        os.remove(simulation.output_file)

    simulation.run(n_steps=10000, timestep=1e-3)

    traj = readdy.Trajectory(simulation.output_file)
    times, counts = traj.read_observable_number_of_particles()
    times = np.array(times) * 1e-3

    stride = 1
    times = times[::stride]
    counts = counts[::stride]

    plt.plot(times, counts[:, 0] + 2. * counts[:, 1], label="A")
    plt.plot(times, counts[:, 2], label="B")
    plt.legend(loc="best")
    plt.show()

    plt.plot(counts[:, 0] + 2. * counts[:, 1], counts[:, 2])
    plt.show()
