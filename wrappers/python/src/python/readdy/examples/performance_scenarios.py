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
Measure the computation time of readdy and its sub-components for pre-defined scenarios.

## Some general ideas for comparing scenarios:

We define the following unit-less quantities (N,L,S,R), which will help to compare different scenarios:
- N ~ number of particles. For reactive systems this is the initial number of particles. Since we do measure for a few timesteps
  only, this is close to the average number of particles.
- L = V / h**3 ~ system size, i.e. fraction of system volume and the largest relevant lengthscale. For a cubic system V**1/3 would
  be the edge length of the box. With repulsing particles of radius r, the largest relevant lengthscale
  would be h=2r. V is the volume _accessible_ to particles. Note that h is the largest distance of two particles
  when an interaction has to be considered. Thus h can also be a reaction distance.
- S = sqrt(2*D*Delta t) / h ~ relative displacement of particles, fraction of distance traveled within
  one timestep (of width Delta t) divided by the largest relevant lengthscale. D is the largest diffusion constant.
  Usually S should be smaller than 1.
- R = k*Delta t ~ reactivity. When k is the largest reaction rate, R should usually be smaller than 1.

Derived from N and L can be the dimensionless density rho = N / L. When there is only a single relevant lengthscale h
or all lengthscales are narrowly distributed around h, the following statements apply:
- when rho = 1, the system is rather dense (every particle is assigned with a cube of volume h**3)
- when rho ~= 2, the system is closely packed (every particles is a sphere of radius h/2)
- when rho < 1, the system is dilute
- when rho >> 1, either lengthscales are widely distributed (then rho is ambiguous) or the system is ill-parametrized

There are different combinations of scenarios where our algorithms can scale very differently:
- Purely reactive system vs. purely collisive systemvs. both reactive and collisive
- Uniformly distributed vs. spatially separated
- Homogeneously sized particles vs. heterogene mixture

If you consider ONE system (fixed N and L) and want to compare different simulation algorithms with
varying timestep Delta t, then the simulation efficiency is of central interest:
- eta = simulated time / computation time ~ simulation efficieny

For example increasing the timestep would increase eta but R and S will also increase. If R or S come into range
of ~1, this hints towards a sampling problem.


## How the following performance scenarios work:

- Set of common variables along which to sample the performance, defined in 'factor_keys'
- Different methods: kernel, reaction_handler are given to scenario as parameters
- Scenario performs the system configuration and makes use of common variable factors, also runs the simulation
- Sampler that samples any scenario across any given variable (in form of factors from the base value)
"""

import numpy as np
import matplotlib.pyplot as plt

import readdy._internal as api

__author__ = "chrisfroe"
__license__ = "LGPL"

# these are varied by the user
factor_keys = ["n_particles", "box_length", "diff", "skin", "reaction_rate"]
# these usually show up on the y-axis
result_keys = ["/", "/integrator", "/neighborList", "/forces", "/reactionScheduler"]
# these usually show up on the x-axis, they are unit-less quantities calculated by the scenario from the factors and hardcoded base values
system_variables_keys = ["n_particles", "system_size", "displacement", "reactivity", "density"]


def get_empty_result_container(n_samples):
    container = dict()
    for key in result_keys:
        container[key] = np.zeros(n_samples)
    return container


def get_identity_factors():
    identity_factors = dict()
    for key in factor_keys:
        identity_factors[key] = 1.
    return identity_factors


def get_empty_system_variables(n_samples):
    system_vars = dict()
    for key in system_variables_keys:
        system_vars[key] = np.zeros(n_samples)
    return system_vars


class PerformanceScenario:
    """
    PerformanceScenario is a thin wrapper for Simulation. Derived classes take
    a dictionary of factors (that scale typical variables like number of particles),
    and configure a certain scenario. The scenario is then run, and performance
    times are set or appended to containers (which usually are dictionaries of lists/arrays).
    """

    def __init__(self, kernel, time_step, integrator, reaction_scheduler):
        self.sim = api.Simulation()
        self.sim.set_kernel(kernel)
        if integrator is not None:
            self.integrator = integrator
        else:
            self.integrator = "EulerBDIntegrator"
        if reaction_scheduler is not None:
            self.reaction_scheduler = reaction_scheduler
        else:
            self.reaction_scheduler = "UncontrolledApproximation"

        self.time_step = time_step
        self.system_vars = dict()
        for key in system_variables_keys:
            self.system_vars[key] = 0.

    def run(self, n_steps, **kwargs):
        conf = self.sim.run_scheme_readdy()
        conf.with_integrator(self.integrator)
        conf.with_reaction_scheduler(self.reaction_scheduler)
        conf.configure_and_run(n_steps, self.time_step)
        if "skin" in kwargs:
            conf.with_skin_size(kwargs["skin"])

    def performance(self):
        return self.sim.performance_root()

    def set_times(self, container, idx):
        """Set results with specific index if container elements are array-like"""
        for key in result_keys:
            container[key][idx] = self.performance()[key].time()

    def set_counts(self, container, idx):
        """Set results with specific index if container elements are array-like"""
        for key in result_keys:
            container[key][idx] = self.performance()[key].count()

    def set_system_vars(self, container, idx):
        """Set system variables with specific index if container elements are array-like"""
        for key in system_variables_keys:
            container[key][idx] = self.system_vars[key]

    @classmethod
    def describe(cls):
        return "PerformanceScenario"


class Collisive(PerformanceScenario):
    """Scenario with uniformly distributed particles that repulse each other"""

    def __init__(self, kernel, factors, time_step=0.01, integrator=None, reaction_scheduler=None):
        super(Collisive, self).__init__(kernel, time_step, integrator, reaction_scheduler)
        box_length = 7. * factors["box_length"]
        self.sim.box_size = api.Vec(box_length, box_length, box_length)
        self.sim.periodic_boundary = [True, True, True]
        diffusion_coeff = 1. * factors["diff"]
        self.sim.register_particle_type("A", diffusion_coeff, 0.5)
        interaction_distance = 1.  # currently due to particle radius 0.5
        self.sim.register_potential_harmonic_repulsion("A", "A", 10.)

        # set dimensionless system quantities
        self.system_vars["n_particles"] = int(100 * factors["n_particles"])
        self.system_vars["system_size"] = box_length ** 3 / interaction_distance ** 3
        self.system_vars["displacement"] = np.sqrt(2. * diffusion_coeff * self.time_step) / interaction_distance
        self.system_vars["reactivity"] = 0.
        self.system_vars["density"] = self.system_vars["n_particles"] / self.system_vars["system_size"]

        for i in range(self.system_vars["n_particles"]):
            pos = np.random.uniform(size=3) * box_length - 0.5 * box_length
            self.sim.add_particle("A", api.Vec(*pos))

    @classmethod
    def describe(cls):
        return "Collisive"


class Reactive(PerformanceScenario):
    """Scenario with three species uniformly distributed and reactions"""

    def __init__(self, kernel, factors, time_step=0.01, integrator=None, reaction_scheduler=None):
        super(Reactive, self).__init__(kernel, time_step, integrator, reaction_scheduler)
        box_length = 20. * factors["box_length"]
        self.sim.box_size = api.Vec(box_length, box_length, box_length)
        self.sim.periodic_boundary = [True, True, True]
        diffusion_coeff = 1. * factors["diff"]
        self.sim.register_particle_type("A", diffusion_coeff, 0.5)
        self.sim.register_particle_type("B", diffusion_coeff, 0.5)
        self.sim.register_particle_type("C", diffusion_coeff, 0.5)
        reaction_radius = 2.
        rate_on = 1e-3 * factors["reaction_rate"]
        rate_off = 5e-5 * factors["reaction_rate"]
        self.sim.register_reaction_fusion("fusion", "A", "B", "C", rate_on, reaction_radius)
        self.sim.register_reaction_fission("fission", "C", "A", "B", rate_off, reaction_radius)
        n_a = int(25 * factors["n_particles"])
        n_c = int(90 * factors["n_particles"])

        # set dimensionless system quantities
        self.system_vars["n_particles"] = 2 * n_a + n_c
        self.system_vars["system_size"] = box_length ** 3 / reaction_radius ** 3
        self.system_vars["displacement"] = np.sqrt(2. * diffusion_coeff * self.time_step) / reaction_radius
        self.system_vars["reactivity"] = rate_on * time_step
        self.system_vars["density"] = self.system_vars["n_particles"] / self.system_vars["system_size"]

        for i in range(self.system_vars["n_particles"]):
            pos = np.random.uniform(size=3) * box_length - 0.5 * box_length
            self.sim.add_particle("A", api.Vec(*pos))
            pos = np.random.uniform(size=3) * box_length - 0.5 * box_length
            self.sim.add_particle("B", api.Vec(*pos))
            pos = np.random.uniform(size=3) * box_length - 0.5 * box_length
            self.sim.add_particle("C", api.Vec(*pos))

    @classmethod
    def describe(cls):
        return "Reactive"


class ReactiveCollosive(PerformanceScenario):
    # todo combine both
    pass


class Sampler:
    """
    Samplers only contain a function sample(). For plotting it is convenient to
    pass this function around together with a describe() method.
    """
    @staticmethod
    def sample(*args, **kwargs):
        pass

    @classmethod
    def describe(cls):
        return "Sampler"


class ConstDensity(Sampler):
    @staticmethod
    def sample(number_factors, n_time_steps=50, scenario_type=Collisive, kernel="SingleCPU"):
        """Helper function to sample a scenario for varying number of particles at constant density

        :param number_factors: array of factors that will be applied to the scenarios, the length of this determines the number of samples/simulations
        :param n_time_steps: how many time steps shall be performed for each sample/simulation
        :param scenario_type: determines which scenario shall be run,
        :param kernel: string determining the compute kernel of the simulation
        :return: performance results (times, counts) and unit-less system variables, these are all dictionaries
        """
        n_samples = len(number_factors)
        times = get_empty_result_container(n_samples)
        counts = get_empty_result_container(n_samples)
        system_vars = get_empty_system_variables(n_samples)
        box_length_factors = np.cbrt(number_factors)
        for i, _ in enumerate(number_factors):
            factors = get_identity_factors()
            factors["n_particles"] = number_factors[i]
            factors["box_length"] = box_length_factors[i]
            scenario = scenario_type(kernel, factors)
            scenario.run(n_time_steps)
            scenario.set_system_vars(system_vars, i)
            scenario.set_times(times, i)
            scenario.set_counts(counts, i)
        return times, counts, system_vars

    @classmethod
    def describe(cls):
        return "ConstDensity"


class ConstVolume(Sampler):
    @staticmethod
    def sample(number_factors, n_time_steps=50, scenario_type=Collisive, kernel="SingleCPU"):
        """Helper function to sample a scenario for varying number of particles at constant volume

        :param number_factors: array of factors that will be applied to the scenarios, the length of this determines the number of samples/simulations
        :param n_time_steps: how many time steps shall be performed for each sample/simulation
        :param scenario_type: determines which scenario shall be run,
        :param kernel: string determining the compute kernel of the simulation
        :return: performance results (times, counts) and unit-less system variables, these are all dictionaries
        """
        n_samples = len(number_factors)
        times = get_empty_result_container(n_samples)
        counts = get_empty_result_container(n_samples)
        system_vars = get_empty_system_variables(n_samples)
        box_length_factors = np.ones_like(number_factors) * 2.
        for i, _ in enumerate(number_factors):
            factors = get_identity_factors()
            factors["n_particles"] = number_factors[i]
            factors["box_length"] = box_length_factors[i]
            scenario = scenario_type(kernel, factors)
            scenario.run(n_time_steps)
            scenario.set_system_vars(system_vars, i)
            scenario.set_times(times, i)
            scenario.set_counts(counts, i)
        return times, counts, system_vars

    @classmethod
    def describe(cls):
        return "ConstVolume"


def plot_times(times, counts, x_axis, total_label=None):
    """Helper function to plot results of a sampling

    :param times: performance results dictionary
    :param counts: performance results dictionary
    :param x_axis: an array, most likely one of system_vars
    :param total_label: label if only total time shall be plotted
    :return:
    """
    if total_label is not None:
        plt.plot(x_axis, times["/"], label=total_label)
        plt.ylabel("total computation time in seconds")
    else:
        plt.plot(x_axis, times["/integrator"] / counts["/integrator"], label="integrator")
        plt.plot(x_axis, times["/neighborList"] / counts["/neighborList"], label="neighbor-list, runs twice")
        plt.plot(x_axis, times["/forces"] / counts["/forces"], label="forces")
        plt.plot(x_axis, times["/reactionScheduler"] / counts["/reactionScheduler"], label="reactions")
        plt.ylabel("computation time per evaluation in seconds")
