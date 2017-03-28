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


from __future__ import print_function
from readdy._internal.readdybinding.api import KernelProvider, Simulation
from readdy._internal.readdybinding.common import Vec

import numpy as np

from readdy.util import platform_utils

import matplotlib.pyplot as plt

rdf = None
centers = None
n_calls = 0
# T = 200000000
T = 20000000


def rdf_callback(pair):
    global rdf, centers, n_calls, T
    if centers is None:
        centers = pair[0]
    if rdf is None:
        rdf = pair[1]
    else:
        rdf += pair[1]
    n_calls += 1
    if n_calls % 10000 == 0:
        print("%s" % (10.*float(n_calls)/float(T)))


if __name__ == '__main__':
    KernelProvider.get().load_from_dir(platform_utils.get_readdy_plugin_dir())
    simulation = Simulation()
    simulation.set_kernel("CPU")

    box_size = Vec(10, 10, 10)
    simulation.kbt = 2
    simulation.periodic_boundary = [True, True, True]
    simulation.box_size = box_size
    simulation.register_particle_type("A", .2, 1.)
    simulation.register_particle_type("B", .2, 1.)
    simulation.register_potential_harmonic_repulsion("A", "B", 10)
    simulation.add_particle("A", Vec(-2.5, 0, 0))
    simulation.add_particle("B", Vec(0, 0, 0))

    simulation.register_observable_radial_distribution(10, np.arange(0, 5, .01), ["A"], ["B"], 1. / (box_size[0] * box_size[1] * box_size[2]), rdf_callback)
    simulation.run(T, 0.02)

    print("n_calls=%s" % n_calls)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(centers, rdf / n_calls)
    plt.show()
