#!/usr/bin/env python
# coding=utf-8

# Copyright © 2017 Computational Molecular Biology Group,
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
Execute performance measurement and generate some plots
"""

import numpy as np
import matplotlib.pyplot as plt

import readdy.examples.performance_scenarios as ps

__author__ = "chrisfroe"
__license__ = "LGPL"

if __name__ == '__main__':
    n_samples = 25
    n_time_steps = 50
    results = {"SingleCPU": None, "CPU": None}
    for kernel in ["SingleCPU", "CPU"]:
        number_factors = np.logspace(0, 3.5, n_samples)
        times, counts, system_vars = ps.sample_n_particles_const_density(number_factors, n_time_steps, scenario_type=ps.Collisive, kernel=kernel)
        results[kernel] = (times, counts, system_vars)
        print("density", system_vars["density"])
        print("reactivity", system_vars["reactivity"])
        print("displacement", system_vars["displacement"])
        ps.plot_times(times, counts, system_vars["n_particles"])
        plt.legend(loc="best")
        plt.title("Repulsing particles at constant density, " + kernel)
        plt.xlabel("number of particles")
        plt.xscale("log")
        plt.yscale("log")
        plt.savefig("repulsion-"+kernel+".pdf")
        plt.clf()

    for kernel in ["SingleCPU", "CPU"]:
        ps.plot_times(results[kernel][0], results[kernel][1], results[kernel][2]["n_particles"], total=kernel)

    plt.title("Repulsing particles at constant density, SCPU vs CPU")
    plt.xlabel("number of particles")
    plt.xscale("log")
    plt.yscale("log")
    plt.legend(loc="best")
    plt.savefig("repulsion-SCPU-vs-CPU.pdf")