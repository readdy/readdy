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
Created on 08.09.17

@author: clonker, chrisfroe
"""
from readdy.api.conf.KernelConfiguration import CPUKernelConfiguration as _CPUKernelConfiguration
from readdy.api.conf.KernelConfiguration import NOOPKernelConfiguration as _NOOPKernelConfiguration
from readdy.api.registry.observables import Observables as _Observables
from readdy._internal.readdybinding.api import Simulation as _Simulation
from readdy._internal.readdybinding.api import Saver as _Saver
from readdy._internal.readdybinding.api import UserDefinedAction as _UserDefinedAction
from readdy.api.utils import vec3_of as _v3_of
from readdy._internal.readdybinding.api import KernelProvider as _KernelProvider

from readdy.api.experimental.action_factory import ActionFactory as _ActionFactory


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
        self._simulation = _Simulation(kernel, context)

        self._evaluate_topology_reactions = True
        self._evaluate_forces = True
        self._evaluate_observables = True
        self._skin = 0.
        self._integrator = "EulerBDIntegrator"
        self._reaction_handler = "Gillespie"
        self._simulation_scheme = "ReaDDyScheme"
        self._output_file = ""

        self.output_file = output_file
        self._observables = _Observables(self)
        # fixme self._checkpoint_saver = None
        self._make_checkpoints = False
        self._checkpoint_format = "checkpoint_{}.h5"
        self._checkpoint_stride = None
        self._checkpoint_outdir = None
        self._checkpoint_max_n_saves = 5

        self.integrator = integrator
        self.reaction_handler = reaction_handler
        self.evaluate_topology_reactions = evaluate_topology_reactions
        self.evaluate_forces = evaluate_forces
        self.evaluate_observables = evaluate_observables
        self.skin = skin
        self._show_progress = True
        self._progress_output_stride = 10

        if kernel == "CPU":
            self._kernel_configuration = _CPUKernelConfiguration()
        else:
            self._kernel_configuration = _NOOPKernelConfiguration()

        self._action_factory = _ActionFactory(self._simulation)

    @property
    def _actions(self):
        return self._action_factory

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
    def progress_output_stride(self) -> int:
        """
        If the progress bar is shown, this value determines how often it is updated. A value of n means that every
        n-th time step the progress is increased by one.
        :return: the progress output stride
        """
        return self._progress_output_stride

    @progress_output_stride.setter
    def progress_output_stride(self, value: int):
        """
        Sets the progress output stride, determining how often the progressbar (if any) is updated. A value of n means
        that every n-th time step the progress is increased by one.
        :param value: the progress output stride
        """
        assert value > 0, "The progress output stride must be a positive value but was {}".format(value)
        self._progress_output_stride = value

    @property
    def show_progress(self) -> bool:
        """
        Returns if a progress bar is shown during simulation. Will only appear if there are more than 100 time steps.
        :return: true if a progress bar should be shown
        """
        return self._show_progress

    @show_progress.setter
    def show_progress(self, value: bool):
        """
        Sets whether to show a progress bar.
        :param value: true if a progress bar should be shown
        """
        self._show_progress = value

    @property
    def evaluate_topology_reactions(self) -> bool:
        """
        Returns whether to evaluate topology reactions during simulation.
        :return: a boolean
        """
        return self._evaluate_topology_reactions

    @evaluate_topology_reactions.setter
    def evaluate_topology_reactions(self, value: bool):
        """
        Sets whether to evaluate topology reactions during simulation.
        :param value: a boolean value
        """
        assert isinstance(value, bool), "the value must be bool but was {}".format(type(value))
        self._evaluate_topology_reactions = value

    @property
    def evaluate_forces(self) -> bool:
        """
        Returns whether to evaluate forces during simulation.
        :return: a boolean
        """
        return self._evaluate_forces

    @evaluate_forces.setter
    def evaluate_forces(self, value: bool):
        """
        Sets whether to evaluate forces during simulation.
        :param value: a boolean value
        """
        assert isinstance(value, bool), "the value must be bool but was {}".format(type(value))
        self._evaluate_forces = value

    @property
    def evaluate_observables(self) -> bool:
        """
        Returns whether observables are evaluated during simulation.
        :return: a boolean
        """
        return self._evaluate_observables

    @evaluate_observables.setter
    def evaluate_observables(self, value: bool):
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
        assert (isinstance(value, str) and value in supported_integrators) or isinstance(value, _UserDefinedAction), \
            "the integrator can only be one of {} or a user defined integrator.".format(",".join(supported_integrators))
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
            * DetailedBalance
        :param value: the reaction handler
        """
        supported_reaction_handlers = ("Gillespie", "UncontrolledApproximation", "DetailedBalance",)
        assert (isinstance(value, str) and value in supported_reaction_handlers) or \
               isinstance(value, _UserDefinedAction), \
            "the reaction handler can only be one of {} or a user defined reaction handler" \
                .format(",".join(supported_reaction_handlers))
        self._reaction_handler = value

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

    def make_checkpoints(self, stride, output_directory, max_n_saves=5):
        """
        Records the system's state (particle positions and topology configuration) every stride steps into the
        trajectory file. This can be used to load particle positions to continue a simulation.

        :param stride: record a checkpoint every `stride` simulation steps
        :param output_directory: directory containing checkpoint files
        :param max_n_saves: only keep `max_n_saves` many checkpoint files, in case of `max_n_saves=0` all files are kept
        """
        import os
        if not os.path.exists(output_directory):
            os.makedirs(output_directory)
        self._checkpoint_outdir = output_directory
        self._checkpoint_max_n_saves = max_n_saves
        # fixme self._checkpoint_saver = _Saver(str(output_directory), max_n_saves, "checkpoint_{}.h5")
        self._checkpoint_stride = stride
        self._make_checkpoints = True

    @staticmethod
    def list_checkpoints(file_name):
        """
        Returns a list of checkpoints in a trajectory file.
        :param file_name: the trajectory file
        """
        from readdy import Trajectory
        traj = Trajectory(file_name)
        return traj.list_checkpoints()

    @staticmethod
    def list_checkpoint_files(directory_name):
        import os
        from glob import glob
        from readdy import Trajectory
        files = glob(os.path.join(directory_name, '*.h5'))
        result = []
        for f in files:
            try:
                traj = Trajectory(f)
                traj.list_checkpoints()
                result.append(f)
            except:
                pass
        return result

    @staticmethod
    def get_latest_checkpoint_file(directory_name):
        """
        Yields the file with the most up-to-date checkpoint in it in terms of the 'step'.
        :param directory_name: the directory in which to look for checkpoints
        :return: the latest checkpoint or raise if directory is empty
        """
        from readdy import Trajectory

        checkpoint_files = Simulation.list_checkpoint_files(directory_name)
        if len(checkpoint_files) == 0:
            raise ValueError("No checkpoints found in {}".format(directory_name))
        latest = None
        latest_step = None
        for filename in checkpoint_files:

            traj = Trajectory(filename)
            ckpt = traj.list_checkpoints()

            # print(f"GetLatestCheckpointFile: File {filename} checkpoints {ckpt}")

            latest_ckpt_step = 0
            for c in ckpt:
                latest_ckpt_step = c['step'] if latest_step is None else (c['step'] if c['step'] > latest_step
                                                                          else latest_step)
            if latest is None:
                latest = (filename, ckpt)
                latest_step = latest_ckpt_step
                # print(f"  -> got new latest {latest}")
            else:
                if latest_ckpt_step > latest_step:
                    latest = (filename, ckpt)
                    latest_step = latest_ckpt_step
                    # print(f"  -> got new latest {latest}")
        return latest[0]

    def load_particles_from_latest_checkpoint(self, directory_name):
        """
        Loads the latest checkpoint in terms of 'step'.
        :param directory_name: the directory in which to look for checkpoint files
        """
        f = self.get_latest_checkpoint_file(directory_name)
        self.load_particles_from_checkpoint(f, n=None)

    def load_particles_from_checkpoint(self, file_name, n=None):
        """
        Adds particles to the simulation as contained in the n-th checkpoint of a trajectory file.
        :param file_name: the trajectory file
        :param n: if n is None, retrieve configuration from latest checkpoint, otherwise use 'n-th' checkpoint, n >= 0
        """
        import numpy as _np
        from readdy import Trajectory
        from readdy.api.trajectory import _CKPT
        from readdy.util.io_utils import get_particle_types
        checkpoints = self.list_checkpoints(file_name)
        if n is None:
            n = len(checkpoints)-1
        else:
            assert n < len(checkpoints), f"n={n} is out of bounds, only have {len(checkpoints)} checkpoints"
        assert n >= 0, f"n must be positive but was {n}"

        # group particle types by flavor (NORMAL, TOPOLOGY)
        ptypes = get_particle_types(filename=file_name)
        normal_types = []
        topology_types = []
        for t in ptypes.values():
            if t['flavor'] == 'TOPOLOGY':
                topology_types.append(t['type_id'])
            if t['flavor'] == 'NORMAL':
                normal_types.append(t['type_id'])

        # load frame into memory
        traj = Trajectory(file_name)
        n_particles_per_frame, positions, types, ids = traj.to_numpy(start=n, stop=n+1, name=_CKPT.POSITIONS_CKPT)

        # add particles with flavor NORMAL
        for normal_type in normal_types:
            tixs = _np.argwhere(types[0] == normal_type)
            pos = positions[0][tixs].squeeze()
            if len(pos) > 0:
                self.add_particles(traj.species_name(normal_type), _np.atleast_2d(pos))

        # add topologies
        time, topology_records = traj.read_observable_topologies(start=n, stop=n+1, data_set_name=_CKPT.TOPOLOGY_CKPT)
        assert len(topology_records) == 1
        for topology in topology_records[0]:
            particle_types = [traj.species_name(types[0, i]) for i in topology.particles]
            pos = _np.atleast_2d(_np.array([positions[0, i] for i in topology.particles]))
            top = self.add_topology(traj.topology_type_name(topology.type), particle_types, pos)
            for e in topology.edges:
                top.graph.add_edge(e[0], e[1])

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
        assert positions.shape[1] == 3, "shape[1] has to be 3 but shape was {}".format(positions.shape)
        if isinstance(particle_types, str):
            particle_types = [particle_types]
        return self._simulation.add_topology(topology_type, particle_types, positions)

    @property
    def current_particles(self):
        """
        Returns a list of currently added particles.
        :return: a list of currently added particles
        """
        return self._simulation.current_particles

    @property
    def current_topologies(self):
        """
        Returns a list of currently registered topologies
        :return: a list of topologies
        """
        return self._simulation.current_topologies

    def run(self, n_steps, timestep, show_summary=True):
        """
        Executes the simulation as configured.

        :param n_steps: number of steps to perform
        :param timestep: the time step to use [time]
        :param show_summary: determines if system and simulation configuration is printed
        """
        import os
        from contextlib import closing
        import readdy._internal.readdybinding.common.io as io
        # from readdy._internal.readdybinding.common.util import ostream_redirect

        self._simulation.context.validate()

        timestep = self._unit_conf.convert(timestep, self.time_unit)
        assert timestep > 0.

        if self.output_file is not None and len(self.output_file) > 0 and os.path.exists(self.output_file):
            raise ValueError("Output file already existed: {}".format(self.output_file))

        self._simulation.set_kernel_config(self.kernel_configuration.to_json())

        loop = self._simulation.create_loop(timestep)
        loop.use_integrator(self.integrator)
        loop.evaluate_forces(self.evaluate_forces)
        loop.evaluate_topology_reactions(self.evaluate_topology_reactions, timestep)
        loop.use_reaction_scheduler(self.reaction_handler)
        loop.evaluate_observables(self.evaluate_observables)
        if self.integrator == "MdgfrdIntegrator":
            loop.neighbor_list_cutoff = max(2. * self._simulation.context.calculate_max_cutoff(), loop.neighbor_list_cutoff)
        if self._skin > 0.:
            loop.neighbor_list_cutoff = loop.neighbor_list_cutoff + self._skin
        if self._make_checkpoints:
            loop.make_checkpoints(self._checkpoint_stride, self._checkpoint_outdir,
                                  self._checkpoint_max_n_saves, self._checkpoint_format)

        write_outfile = self.output_file is not None and len(self.output_file) > 0

        import contextlib

        @contextlib.contextmanager
        def nullcontext():
            yield

        with closing(io.File.create(self.output_file)) if write_outfile else nullcontext() as f:
            if write_outfile:
                for name, chunk_size, handle in self._observables._observable_handles:
                    handle.enable_write_to_file(f, name, chunk_size)
                loop.write_config_to_file(f)
            if show_summary:
                print(self._simulation.context.describe())
                print(loop.describe())
            if self.show_progress:
                import tqdm
                loop.progress_output_stride = self._progress_output_stride
                with tqdm.tqdm(total=n_steps // self._progress_output_stride) as progress:
                    loop.progress_callback = lambda t: progress.update(1)
                    loop.run(n_steps)
                    loop.progress_callback = None
            else:

                loop.run(n_steps)

    def _run_custom_loop(self, custom_loop_function, show_summary=True):
        """
        Executes the simulation loop provided by argument, additionally takes care of preparing output file.

        :param custom_loop_function: function to be executed, should contain calls to actions
        :param show_summary: determines if system configuration is printed
        """
        import os
        from contextlib import closing
        import readdy._internal.readdybinding.common.io as io

        self._simulation.context.validate()

        if self.output_file is not None and len(self.output_file) > 0 and os.path.exists(self.output_file):
            raise ValueError("Output file already existed: {}".format(self.output_file))

        self._simulation.set_kernel_config(self.kernel_configuration.to_json())

        # todo consider exposing enable_write_to_file as method of simulation (cpp and py),
        # todo then we can delete this whole run_custom_loop method
        # the actions created by the loop internally are not used
        # the loop object here is merely used as convenient access to configuration
        # e.g. enable_write_to_file, write_config_to_file, run_initialize
        loop = self._simulation.create_loop(42.)

        write_outfile = self.output_file is not None and len(self.output_file) > 0

        import contextlib

        @contextlib.contextmanager
        def nullcontext():
            yield

        with closing(io.File.create(self.output_file)) if write_outfile else nullcontext() as f:
            if write_outfile:
                for name, chunk_size, handle in self._observables._observable_handles:
                    handle.enable_write_to_file(f, name, chunk_size)
                loop.write_config_to_file(f)
                loop.run_initialize()  # writes the config to file here
            if show_summary:
                print(self._simulation.context.describe())
            custom_loop_function()
