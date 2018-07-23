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
Created on 21.06.17

@author: clonker
"""

from __future__ import print_function

from contextlib import closing

import numpy as np
import readdy._internal.readdybinding.api.top as top
import readdy._internal.readdybinding.common.io as io
import readdy._internal.readdybinding.common as common
from readdy._internal.readdybinding.api import KernelProvider
from readdy._internal.readdybinding.api import ParticleTypeFlavor
from readdy._internal.readdybinding.api import Simulation
from readdy.util import platform_utils

class ChainDecay(object):

    def __init__(self, kernel, time_step):
        self.kernel_provider = KernelProvider.get()
        self.kernel_provider.load_from_dir(platform_utils.get_readdy_plugin_dir())
        self.kernel = kernel
        self.time_step = time_step

    def _get_split_reaction(self):
        def reaction_function(topology):
            recipe = top.Recipe(topology)
            if topology.get_n_particles() > 1:
                edge = np.random.randint(0, topology.get_n_particles() - 1)
                recipe.remove_edge(edge, edge + 1)
            return recipe

        def rate_function(topology):
            if topology.get_n_particles() > 1:
                return float(topology.get_n_particles()) / 5.
            else:
                return .0

        fun1 = top.ReactionFunction(reaction_function)
        fun2 = top.RateFunction(rate_function)

        reaction = top.TopologyReaction(fun1, fun2)
        reaction.roll_back_if_invalid()
        reaction.create_child_topologies_after_reaction()
        return reaction

    def _get_decay_reaction(self, typeidb):
        def reaction_function(topology):
            recipe = top.Recipe(topology)
            if topology.get_n_particles() == 1:
                recipe.change_particle_type(0, typeidb)
            return recipe

        def rate_function(topology):
            return 1./self.time_step if topology.get_n_particles() == 1 else 0

        fun1, fun2 = top.ReactionFunction(reaction_function), top.RateFunction(rate_function)
        reaction = top.TopologyReaction(fun1, fun2)
        reaction.raise_if_invalid()
        reaction.create_child_topologies_after_reaction()
        return reaction

    def run(self, time_steps, out_file):
        sim = Simulation()
        sim.set_kernel(self.kernel)
        sim.box_size = common.Vec(60, 20, 20)
        sim.periodic_boundary = [True, True, True]

        typeid_b = sim.register_particle_type("B", 1.0, 1.0, ParticleTypeFlavor.NORMAL)
        sim.register_particle_type("Topology A", .5, .5, ParticleTypeFlavor.TOPOLOGY)

        sim.register_potential_harmonic_repulsion("Topology A", "Topology A", 10)
        sim.register_potential_harmonic_repulsion("Topology A", "B", 10)
        sim.register_potential_harmonic_repulsion("B", "B", 10)

        sim.configure_topology_bond_potential("Topology A", "Topology A", 10, 1.)
        sim.configure_topology_angle_potential("Topology A", "Topology A", "Topology A", 10, np.pi)
        # sim.configure_topology_dihedral_potential("Topology A", "Topology A", "Topology A", "Topology A", 1, 1, -np.pi)

        n_elements = 50.
        particles = [sim.create_topology_particle("Topology A", common.Vec(-25. + i, 0, 0))
                     for i in range(int(n_elements))]
        topology = sim.add_topology(particles)

        for i in range(int(n_elements - 1)):
            topology.get_graph().add_edge(i, i + 1)

        topology.add_reaction(self._get_decay_reaction(typeid_b))
        topology.add_reaction(self._get_split_reaction())

        traj_handle = sim.register_observable_flat_trajectory(1)
        with closing(io.File(out_file, io.FileAction.CREATE, io.FileFlag.OVERWRITE)) as f:
            traj_handle.enable_write_to_file(f, u"", 50)
            sim.run_scheme_readdy(True)\
                .evaluate_topology_reactions()\
                .write_config_to_file(f)\
                .configure_and_run(time_steps, self.time_step)
        print("currently %s topologies" % len(sim.current_topologies()))


if __name__ == '__main__':
    sim = ChainDecay("SingleCPU", .001)
    sim.run(10000, "out.h5")