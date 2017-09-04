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
Measure the computation time of readdy and its
sub-components for pre-defined scenarios.

We define the following unitless quantities (N,L,S,R), which will help to compare different scenarios:
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

 - Derived from N and L can be the dimensionless density rho = N / L. When there is only a single relevant lengthscale h
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
"""

import numpy as np
import readdy._internal.readdybinding.common.perf as perf
import readdy._internal.readdybinding.api as api

__author__ = "chrisfroe"

# we need
# - set of common variables along which to sample the performance:
#   scalar "skin", int "n_particles", scalar "box_length",
#   scalar "diffusion_coeff", scalar "reaction_rate"
# - different methods: string "kernel", string "reaction_handler" are given to scenario as parameters
# - Scenario that has the system configuration and makes use of common variable factors with respect to base values,
#   also runs the simulation with
# - Sampler that samples any scenario across any given variable (in form of factors from the base value)
#   this need not be an object
#
# usage should be done outside of this file such that data can be readily plotted

factor_keys = ["n_particles", "box_length", "diff", "skin"]


class PerformanceScenario:
    def __init__(self, kernel, integrator, reaction_scheduler):
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

    def run(self, n_steps, time_step):
        conf = self.sim.run_scheme_readdy()
        conf.with_integrator(self.integrator)
        conf.with_reaction_scheduler(self.reaction_scheduler)
        conf.configure_and_run(n_steps, time_step)


class Collisive(PerformanceScenario):
    def __init__(self, kernel, factors=None, integrator=None, reaction_scheduler=None):
        super(Collisive, self).__init__(kernel, integrator, reaction_scheduler)
        box_length = 20. * factors["box_length"]
        self.sim.box_size = api.Vec(box_length, box_length, box_length)
        self.sim.periodic_boundary = [True, True, True]
        self.sim.register_particle_type("A", 1. * factors["diff"], 1.)
        for i in range(100 * factors["n_particles"]):
            pos = np.random.uniform(size=3) * box_length - 0.5 * box_length
            self.sim.add_particle("A", api.Vec(*pos))
