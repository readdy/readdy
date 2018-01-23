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
from readdy.util.progress import SimulationProgress as _SimulationProgress
from readdy._internal.readdybinding.api import KernelProvider as _KernelProvider

class Simulation(object):

    def __init__(self, kernel, context, unit_config, output_file="", integrator="EulerBDIntegrator",
                 reaction_handler="Gillespie", evaluate_topology_reactions=True, evaluate_forces=True,
                 evaluate_observables=True, skin=0):
        """
        Creates a new simulation object
        :param kernel: the kernel to use
        :param context: the parent low level context object
        :param unit_config: the unit configuration
        :param output_file: the output file
        :param integrator: the integrator
        :param reaction_handler: the reaction handler
        :param evaluate_topology_reactions: whether to evaluate topology reactions
        :param evaluate_forces: whether to evaluate forces
        :param evaluate_observables: whether to evaluate observables
        :param skin: the skin size for neighbor lists
        """
        available_kernels = _KernelProvider.get().available_kernels()
        if kernel not in available_kernels:
            raise ValueError("The selected kernel was \"{}\" but only {} "
                             "are available.".format(kernel, ", ".join(['"{}"'.format(x) for x in available_kernels])))
        self._unit_conf = unit_config
        self._kernel = kernel
        self._simulation = _Simulation()
        self._simulation.set_kernel(kernel)
        self._simulation.context = context

        self.output_file = output_file
        self._observables = _Observables(self)

        self.integrator = integrator
        self.reaction_handler = reaction_handler
        self.evaluate_topology_reactions = evaluate_topology_reactions
        self.evaluate_forces = evaluate_forces
        self.evaluate_observables = evaluate_observables
        self.skin = skin
        self.simulation_scheme = "ReaDDyScheme"
        self.show_progress = True

        self._progress = None

        if kernel == "CPU":
            self._kernel_configuration = _CPUKernelConfiguration()
        else:
            self._kernel_configuration = _NOOPKernelConfiguration()

    @property
    def units(self):
        return self._unit_conf.reg

    @property
    def length_unit(self):
        return self._unit_conf.length_unit

    @property
    def energy_unit(self):
        return self._unit_conf.energy_unit

    @property
    def time_unit(self):
        return self._unit_conf.time_unit

    @property
    def show_progress(self):
        """
        Returns if a progress bar is shown during simulation. Will only appear if there are more than 100 time steps.
        :return: true if a progress bar should be shown
        """
        return self._show_progress

    @show_progress.setter
    def show_progress(self, value):
        """
        Sets whether to show a progress bar.
        :param value: true if a progress bar should be shown
        """
        self._show_progress = value

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
        :return: the skin size [length]
        """
        return self._skin * self.length_unit

    @skin.setter
    def skin(self, value):
        """
        Sets the skin size to be used in neighbor lists. If the neighbor list is not a Verlet list / adaptive list,
        it will increase the box size in the corresponding cell-linked list.
        :param value: the new skin size [length]
        """
        value = self._unit_conf.convert(value, self.length_unit)
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
        if value in ("ReaDDyScheme", "AdvancedScheme"):
            self._simulation_scheme = value
        else:
            raise ValueError("Simulation scheme value can only be one of \"ReaDDyScheme\" and \"AdvancedScheme\".")

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

    def record_trajectory(self, stride=1, name="", chunk_size=1000):
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
        :param position: the position (ndarray or tuple or list of length 3) [length]
        """
        position = self._unit_conf.convert(position, self.length_unit)
        self._simulation.add_particle(type, _v3_of(position))

    def add_particles(self, type, positions):
        """
        Adds particles of a certain type to the simulation box.

        :param type: type of the particles
        :param positions: (N, 3)-shaped nd-array of positions [length]
        """
        assert positions.shape[1] == 3, "shape[1] has to be 3 but was {}".format(positions.shape[1])
        positions = self._unit_conf.convert(positions, self.length_unit)
        self._simulation.add_particles(type, positions)

    def add_topology(self, topology_type, particle_types, positions):
        """
        Creates and returns a topology by adding `N` topology particles of given type(s) to the simulation.
        Afterwards the connectivity needs to be set.

        :param topology_type: the topology type
        :param particle_types: either a list of types of length `N` or a single string which is then applied as type
                               for all given positions
        :param positions: (N, 3)-shaped nd-array of positions [length]
        :return: the topology object
        """
        positions = self._unit_conf.convert(positions, self.length_unit)
        assert positions.shape[1] == 3, "shape[1] has to be 3 but was {}".format(positions.shape[0])
        if isinstance(particle_types, str):
            particle_types = [particle_types]
        return self._simulation.add_topology(topology_type, particle_types, positions)

    def run(self, n_steps, timestep, show_system=True):
        """
        Executes the simulation as configured.

        :param n_steps: number of steps to perform
        :param timestep: the time step to use [time]
        :param show_system: determines if system configuration is printed
        """
        import os
        from contextlib import closing
        import readdy._internal.readdybinding.common.io as io

        if show_system:
            print(self._simulation.context.describe())

        timestep = self._unit_conf.convert(timestep, self.time_unit)
        assert timestep > 0.

        if self.output_file is not None and len(self.output_file) > 0 and os.path.exists(self.output_file):
            raise ValueError("Output file already existed: {}".format(self.output_file))

        self._simulation.set_kernel_config(self.kernel_configuration.to_json())

        if self.simulation_scheme == 'ReaDDyScheme':
            conf = self._simulation.run_scheme_readdy(False)
        elif self.simulation_scheme == 'AdvancedScheme':
            conf = self._simulation.run_scheme_advanced(False).include_compartments()
        else:
            raise ValueError("Invalid simulation scheme type: {}".format(self.simulation_scheme))

        conf = conf.with_integrator(self.integrator) \
            .include_forces(self.evaluate_forces) \
            .evaluate_topology_reactions(self.evaluate_topology_reactions) \
            .with_reaction_scheduler(self.reaction_handler) \
            .with_skin_size(self._skin) \
            .evaluate_observables(self.evaluate_observables)

        if self.output_file is not None and len(self.output_file) > 0:
            with closing(io.File.create(self.output_file)) as f:
                for name, chunk_size, handle in self._observables._observable_handles:
                    handle.enable_write_to_file(f, name, chunk_size)
                scheme = conf.write_config_to_file(f).configure(timestep)
                if self.show_progress:
                    self._progress = _SimulationProgress(n_steps // 10)
                    scheme.set_progress_callback(self._progress.callback)
                    scheme.set_progress_output_stride(10)
                scheme.run(n_steps)
        else:
            scheme = conf.configure(timestep)
            if self.show_progress:
                self._progress = _SimulationProgress(n_steps // 10)
                scheme.set_progress_callback(self._progress.callback)
                scheme.set_progress_output_stride(10)
            scheme.run(n_steps)
        if self.show_progress:
            self._progress.finish()
