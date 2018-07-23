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

"""
Execute performance measurement and generate some plots
"""

import numpy as np
import matplotlib.pyplot as plt

import readdy.examples.performance_scenarios as ps

__author__ = "chrisfroe"
__license__ = "LGPL"

data_dir = "./"
#  data_dir = "./performance_results/"


def sample_and_plot_scenario(sampler=ps.Sampler, scenario_type=ps.PerformanceScenario, kernel="SingleCPU"):
    n_samples = 5
    n_time_steps = 100
    number_factors = np.logspace(0, 3, n_samples)
    times, counts, system_vars = sampler.sample(number_factors, n_time_steps, scenario_type=scenario_type, kernel=kernel)
    result = (times, counts, system_vars)
    print("density", system_vars["density"])
    print("reactivity", system_vars["reactivity"])
    print("displacement", system_vars["displacement"])
    ps.plot_times(times, counts, system_vars["n_particles"])
    factor = 1e-6
    plt.plot(system_vars["n_particles"], np.fromiter(map(lambda x: factor * x, system_vars["n_particles"]), dtype=np.float), "--", label="linear")
    plt.legend(loc="best")
    plt.title(scenario_type.describe() + sampler.describe() + kernel)
    plt.xlabel("number of particles")
    plt.xscale("log")
    plt.yscale("log")
    plt.savefig(data_dir + scenario_type.describe() + kernel + sampler.describe() + ".pdf")
    plt.clf()
    return result, kernel + scenario_type.describe() + sampler.describe()


if __name__ == '__main__':
    results_cd = [sample_and_plot_scenario(ps.ConstDensity, ps.Collisive, "SingleCPU"),
                  sample_and_plot_scenario(ps.ConstDensity, ps.Collisive, "CPU"),
                  sample_and_plot_scenario(ps.ConstDensity, ps.Reactive, "SingleCPU"),
                  sample_and_plot_scenario(ps.ConstDensity, ps.Reactive, "CPU")]
    results_cv = [sample_and_plot_scenario(ps.ConstVolume, ps.Collisive, "SingleCPU"),
                  sample_and_plot_scenario(ps.ConstVolume, ps.Collisive, "CPU"),
                  sample_and_plot_scenario(ps.ConstVolume, ps.Reactive, "SingleCPU"),
                  sample_and_plot_scenario(ps.ConstVolume, ps.Reactive, "CPU")]

    for result, label in results_cd:
        ps.plot_times(result[0], result[1], result[2]["n_particles"], total_label=label)
    plt.title("SCPU vs CPU")
    plt.xlabel("number of particles")
    plt.xscale("log")
    plt.yscale("log")
    plt.legend(loc="best")
    plt.savefig(data_dir + "SCPU-vs-CPU-ConstDensity.pdf")
    plt.clf()

    for result, label in results_cv:
        ps.plot_times(result[0], result[1], result[2]["n_particles"], total_label=label)
    plt.title("SCPU vs CPU")
    plt.xlabel("number of particles")
    plt.xscale("log")
    plt.yscale("log")
    plt.legend(loc="best")
    plt.savefig(data_dir + "SCPU-vs-CPU-ConstVolume.pdf")
    plt.clf()
