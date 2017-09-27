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
Created on 08.09.17

@author: clonker
"""
from readdy.api.conf.KernelConfiguration import CPUKernelConfiguration as _CPUKernelConfiguration
from readdy.api.conf.KernelConfiguration import NOOPKernelConfiguration as _NOOPKernelConfiguration
from readdy.api.registry.observables import Observables as _Observables
from readdy._internal.readdybinding.api import Simulation as _Simulation
from readdy.api.utils import vec3_of as _v3_of

class Simulation(object):

    def __init__(self, kernel, context, output_file="", integrator="EulerBDIntegrator", reaction_handler="Gillespie",
                 evaluate_topology_reactions=True, evaluate_forces=True, evaluate_observables=True, skin=0):
        """
        Creates a new simulation object
        :param kernel: the kernel to use
        :param context: the parent low level context object
        :param output_file: the output file
        :param integrator: the integrator
        :param reaction_handler: the reaction handler
        :param evaluate_topology_reactions: whether to evaluate topology reactions
        :param evaluate_forces: whether to evaluate forces
        :param evaluate_observables: whether to evaluate observables
        :param skin: the skin size for neighbor lists
        """
        self._kernel = kernel
        self._simulation = _Simulation()
        self._simulation.set_kernel(kernel)
        self._simulation.context = context

        self._output_file = output_file
        self._observables = _Observables(self)
        self._integrator = integrator
        self._reaction_handler = reaction_handler
        self._evaluate_topology_reactions = evaluate_topology_reactions
        self._evaluate_forces = evaluate_forces
        self._evaluate_observables = evaluate_observables
        self._skin = skin
        self._simulation_scheme = "ReaDDyScheme"

        if kernel == "CPU":
            self._kernel_configuration = _CPUKernelConfiguration()
        else:
            self._kernel_configuration = _NOOPKernelConfiguration()

    @property
    def evaluate_topology_reactions(self):
        """
        Returns whether to evaluate topology reactions during simulation.
        :return: a boolean
        """
        return self._evaluate_topology_reactions

    @evaluate_topology_reactions.setter
    def evaluate_topology_reactions(self, value):
        """
        Sets whether to evaluate topology reactions during simulation.
        :param value: a boolean value
        """
        assert isinstance(value, bool), "the value must be bool but was {}".format(type(value))
        self._evaluate_topology_reactions = value

    @property
    def evaluate_forces(self):
        """
        Returns whether to evaluate forces during simulation.
        :return: a boolean
        """
        return self._evaluate_forces

    @evaluate_forces.setter
    def evaluate_forces(self, value):
        """
        Sets whether to evaluate forces during simulation.
        :param value: a boolean value
        """
        assert isinstance(value, bool), "the value must be bool but was {}".format(type(value))
        self._evaluate_forces = value

    @property
    def evaluate_observables(self):
        """
        Returns whether observables are evaluated during simulation.
        :return: a boolean
        """
        return self._evaluate_observables

    @evaluate_observables.setter
    def evaluate_observables(self, value):
        """
        Sets whether to evaluate observables during simulation.
        :param value: a boolean value
        """
        self._evaluate_observables = value

    @property
    def skin(self):
        """
        Returns the skin size to be used in neighbor lists. If the neighbor list is not a Verlet list / adaptive list,
        it will increase the box size in the corresponding cell-linked list.
        :return: the skin size
        """
        return self._skin

    @skin.setter
    def skin(self, value):
        """
        Sets the skin size to be used in neighbor lists. If the neighbor list is not a Verlet list / adaptive list,
        it will increase the box size in the corresponding cell-linked list.
        :param value: the new skin size
        """
        assert value >= 0, "the value has to be non-negative"
        self._skin = value

    @property
    def integrator(self):
        """
        Returns the selected integrator.

        :return: the integrator
        """
        return self._integrator

    @integrator.setter
    def integrator(self, value):
        """
        Sets the integrator. Currently supported:
            * EulerBDIntegrator

        :param value: the integrator
        """
        supported_integrators = ("EulerBDIntegrator",)
        assert isinstance(value, str) and value in supported_integrators, \
            "the integrator can only be one of {}".format(",".join(supported_integrators))
        self._integrator = value

    @property
    def reaction_handler(self):
        """
        Returns the selected reaction handler.

        :return: the reaction handler
        """
        return self._reaction_handler

    @reaction_handler.setter
    def reaction_handler(self, value):
        """
        Sets the reaction handler. Currently supported:
            * Gillespie
            * UncontrolledApproximation
        :param value: the reaction handler
        """
        supported_reaction_handlers = ("Gillespie", "UncontrolledApproximation",)
        assert isinstance(value, str) and value in supported_reaction_handlers, \
            "the reaction handler can only be one of {}".format(",".join(supported_reaction_handlers))
        self._reaction_handler = value


    @property
    def simulation_scheme(self):
        """
        Returns the currently selected simulation scheme.

        :return: the scheme
        """
        return self._simulation_scheme

    @simulation_scheme.setter
    def simulation_scheme(self, value):
        """
        Set the simulation scheme to use in the run loop. Thus far there is the `DefaultScheme` and the
        `AdvancedScheme`.

        :param value: the scheme
        """
        if value in ("DefaultScheme", "AdvancedScheme"):
            self._simulation_scheme = value
        else:
            raise ValueError("Simulation scheme value can only be one of \"DefaultScheme\" and \"AdvancedScheme\".")

    @property
    def kernel_configuration(self):
        """
        Returns the kernel configuration. If the CPU kernel was selected, this gives ability to select the number
        of threads and neighbor list implementations.

        :return: the kernel configuration
        """
        return self._kernel_configuration

    @property
    def output_file(self):
        """
        Returns the currently set output file

        :return: the output file
        """
        return self._output_file

    @output_file.setter
    def output_file(self, value):
        """
        Sets the output file.

        :param value: the output file
        """
        assert isinstance(value, str), "output file can only be a string"
        self._output_file = value

    @property
    def observe(self):
        """
        Gives the ability to observe things during simulation time.

        :return: a reference to the observable registry
        """
        return self._observables

    def record_trajectory(self, stride=1000, name="", chunk_size=1000):
        """
        Record trajectory into file if file name is given. The trajectory consists out of two data sets, one contains
        all positions contiguously and the other one is two dimensional and contains begin and end indices for each
        time step.

        :param stride: skip `stride` time steps before evaluating the observable again
        :param name: the name under which the trajectory can be found
        :param chunk_size: the chunk size with which it is stored
        """
        handle = self._simulation.register_observable_flat_trajectory(stride)
        self._observables._observable_handles.append((name, chunk_size, handle))

    def add_particle(self, type, position):
        """
        Adds a particle of a certain type to a certain position in the simulation box.

        :param type: the type
        :param position: the position (ndarray or tuple or list of length 3)
        """
        self._simulation.add_particle(type, _v3_of(position))

    def add_particles(self, type, positions):
        """

        :param type:
        :param positions:
        """
        assert positions.shape[0] == 3, "shape[0] has to be 3 but was {}".format(positions.shape[0])


    def run(self, n_steps, timestep):
        """
        Executes the simulation as configured.

        :param n_steps: number of steps to perform
        :param timestep: the time step to use
        """
        from contextlib import closing
        import readdy._internal.readdybinding.common.io as io
        if self.simulation_scheme == 'ReaDDyScheme':
            if self.output_file is not None and len(self.output_file) > 0:
                with closing(io.File.create(self.output_file)) as f:
                    for name, chunk_size, handle in self._observables._observable_handles:
                        handle.enable_write_to_file(f, name, chunk_size)
                    self._simulation.run_scheme_readdy(False) \
                        .write_config_to_file(f) \
                        .with_integrator(self.integrator) \
                        .include_forces(self.evaluate_forces) \
                        .evaluate_topology_reactions(self.evaluate_topology_reactions) \
                        .with_reaction_scheduler(self.reaction_handler) \
                        .evaluate_observables(self.evaluate_observables) \
                        .configure(timestep).run(n_steps)
            else:
                self._simulation.run_scheme_readdy(False) \
                    .with_integrator(self.integrator) \
                    .include_forces(self.evaluate_forces) \
                    .evaluate_topology_reactions(self.evaluate_topology_reactions) \
                    .with_reaction_scheduler(self.reaction_handler) \
                    .evaluate_observables(self.evaluate_observables) \
                    .configure(timestep).run(n_steps)
        elif self.simulation_scheme == 'AdvancedScheme':
            if self.output_file is not None and len(self.output_file) > 0:
                with closing(io.File.create(self.output_file)) as f:
                    for name, chunk_size, handle in self._observables._observable_handles:
                        handle.enable_write_to_file(f, name, chunk_size)
                    self._simulation.run_scheme_advanced(False) \
                        .write_config_to_file(f) \
                        .with_integrator(self.integrator) \
                        .include_compartments() \
                        .include_forces(self.evaluate_forces) \
                        .evaluate_topology_reactions(self.evaluate_topology_reactions) \
                        .with_reaction_scheduler(self.reaction_handler) \
                        .evaluate_observables(self.evaluate_observables) \
                        .configure(timestep).run(n_steps)
            else:
                self._simulation.run_scheme_advanced(False) \
                    .with_integrator(self.integrator) \
                    .include_forces(self.evaluate_forces) \
                    .evaluate_topology_reactions(self.evaluate_topology_reactions) \
                    .include_compartments() \
                    .with_reaction_scheduler(self.reaction_handler) \
                    .evaluate_observables(self.evaluate_observables) \
                    .configure(timestep).run(n_steps)
        else:
            raise ValueError("Invalid simulation scheme type: {}".format(self.simulation_scheme))
