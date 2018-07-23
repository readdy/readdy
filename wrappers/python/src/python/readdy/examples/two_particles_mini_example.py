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


from itertools import combinations, product

from readdy._internal.readdybinding.api import KernelProvider, Simulation
from readdy._internal.readdybinding.common import Vec

from readdy.util import platform_utils
import numpy as np
import matplotlib

# matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class TwoParticlesMiniExample(object):
    def __init__(self):
        KernelProvider.get().load_from_dir(platform_utils.get_readdy_plugin_dir())
        self.simulation = Simulation()
        self.simulation.set_kernel("CPU")

        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111, projection='3d')
        plt.ioff()
        self.fig.show()
        self.prev_pos = {}
        self.current_plot = None

        self.T = 4000000

    def ppos_callback(self, pos):
        plt.cla()
        self.ax.set_xlim([-1, 1])
        self.ax.set_ylim([-1, 1])
        self.ax.set_zlim([-1, 1])
        r = [-.75, .75]
        for s, e in combinations(np.array(list(product(r, r, r))), 2):
            if np.sum(np.abs(s - e)) == r[1] - r[0]:
                self.ax.plot3D(*zip(s, e), color="b")
        pA = self.simulation.get_particle_positions("A")
        pB = self.simulation.get_particle_positions("B")
        if len(pA) == 1 and len(pB) == 1:
            A = pA[0]; B = pB[0]
            self.ax.scatter([A[0]], [A[1]], [A[2]], color="g", s=100)
            self.ax.scatter([B[0]], [B[1]], [B[2]], color="r", s=100)
            self.ax.plot3D([A[0], B[0]], [A[1], B[1]], [A[2], B[2]], color="r")
        pC = self.simulation.get_particle_positions("C")
        if len(pC) == 1:
            C = pC[0]
            self.ax.scatter([C[0]], [C[1]], [C[2]], color="b", s=100)
        plt.pause(.001)

    def start(self):
        box_size = Vec(2.0, 2.0, 2.0)
        depth = 2.
        desired_dist = .25
        force_constant = 4 * depth / (desired_dist * desired_dist)
        no_interaction_dist = 1.5
        self.simulation.kbt = 0.01
        self.simulation.periodic_boundary = [False, False, False]
        self.simulation.box_size = box_size
        self.simulation.register_particle_type("A", .1, .1)
        self.simulation.register_particle_type("B", .01, .1)
        self.simulation.register_particle_type("C", .5, .1)
        self.simulation.register_potential_piecewise_weak_interaction("A", "B", force_constant, desired_dist, depth,
                                                                          no_interaction_dist)  # (force constant, desired dist, depth, no interaction dist)
        self.simulation.register_reaction_fusion("fusion", "A", "B", "C", 1000., .3, .5, .5)
        self.simulation.register_reaction_fission("fission", "C", "A", "B", 1000., .25, .5, .5)
        self.simulation.register_potential_box("A", 100., Vec(-.75, -.75, -.75), Vec(1.5, 1.5, 1.5), False)
        self.simulation.register_potential_box("B", 100., Vec(-.75, -.75, -.75), Vec(1.5, 1.5, 1.5), False)
        self.simulation.register_potential_box("C", 100., Vec(-.75, -.75, -.75), Vec(1.5, 1.5, 1.5), False)
        self.simulation.add_particle("A", Vec(-.0, -.0, -.0))
        self.simulation.add_particle("B", Vec(0.1, 0.1, 0.1))
        self.simulation.register_observable_particle_positions(1, [], self.ppos_callback)

        self.simulation.run(self.T, .0001)


if __name__ == '__main__':
    TwoParticlesMiniExample().start()
