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

import numpy as np
import matplotlib.pyplot as plt
import contextlib as cl
import time
import os
import random

import readdy._internal as api
import readdy._internal.readdybinding.api.top as top
from readdy.util import topology_utils
from readdy.util import io_utils

api.set_logging_level("info", python_console_out=False)


# import readdyviewer as rv

class TopologyFusions(object):

    def run(self):
        # the simulation object
        sim = api.Simulation()
        sim.set_kernel("SingleCPU")

        # no periodic boundary
        sim.periodic_boundary = [True, True, True]
        # room temperature
        sim.kbt = 2.437
        # simulation box size
        box_edge_length = 20.
        sim.box_size = api.Vec(box_edge_length, box_edge_length, box_edge_length)

        tid_middle = sim.register_particle_type("middle", .1, .5, api.ParticleTypeFlavor.TOPOLOGY)
        tid_end = sim.register_particle_type("end", .1, .5, api.ParticleTypeFlavor.TOPOLOGY)
        tid_a = sim.register_particle_type("A", 1.0, .5)

        sim.register_topology_type("TT")

        sim.configure_topology_bond_potential("middle", "middle", 30., 1.)
        sim.configure_topology_bond_potential("middle", "end", 30., 1.)
        sim.configure_topology_bond_potential("end", "end", 30., 1.)
        sim.register_potential_harmonic_repulsion("middle", "middle", force_constant=30.)
        sim.register_potential_harmonic_repulsion("middle", "end", force_constant=30.)
        sim.register_potential_harmonic_repulsion("end", "end", force_constant=30.)
        sim.register_potential_harmonic_repulsion("A", "A", force_constant=30.)
        sim.register_potential_harmonic_repulsion("A", "middle", force_constant=30.)
        sim.register_potential_harmonic_repulsion("A", "end", force_constant=30.)

        chain = [sim.create_topology_particle("end", api.Vec(-1.,0.,0.)),
                 sim.create_topology_particle("middle", api.Vec(0.,0.,0.)),
                 sim.create_topology_particle("end", api.Vec(0.,0.,1.))]
        topology = sim.add_topology("TT", chain)
        topology.get_graph().add_edge(0, 1)
        topology.get_graph().add_edge(1, 2)

        for i in range(10):
            offsetx = 10 * np.random.random() - 5.
            offsety = 10 * np.random.random() - 5.
            chain = [sim.create_topology_particle("end", api.Vec(-1.,offsetx,offsety)),
                     sim.create_topology_particle("middle", api.Vec(0.,offsetx,offsety)),
                     sim.create_topology_particle("end", api.Vec(1.,offsetx,offsety))]
            topology = sim.add_topology("TT", chain)
            topology.get_graph().add_edge(0, 1)
            topology.get_graph().add_edge(1, 2)

        sim.register_spatial_topology_reaction("Merge: TT(end) + TT(end)->TT(middle--middle)", 10., .6)

        # define observables and run
        sim.register_observable_n_particles(100, ["A", "middle", "end"], lambda x: print("currently {} particles".format(x)))

        sim.run_scheme_readdy() \
            .evaluate_topology_reactions() \
            .with_reaction_scheduler("UncontrolledApproximation") \
            .with_skin_size(5.) \
            .configure_and_run(50000, 0.005)

if __name__ == '__main__':
    tf = TopologyFusions()
    tf.run()