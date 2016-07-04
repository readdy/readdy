from __future__ import print_function
from readdy._internal.simulation import KernelProvider, Simulation, Vec

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
    kernel_provider = KernelProvider.get()
    kernel_provider.load_from_dir(platform_utils.get_readdy_plugin_dir())
    simulation = Simulation()
    simulation.setKernel("SingleCPU")

    box_size = Vec(10, 10, 10)
    simulation.kbt = 2
    simulation.periodic_boundary = [True, True, True]
    simulation.box_size = box_size
    simulation.registerParticleType("A", .2, 1.)
    simulation.registerParticleType("B", .2, 1.)
    simulation.registerHarmonicRepulsionPotential("A", "B", 10)
    simulation.addParticle("A", Vec(-2.5, 0, 0))
    simulation.addParticle("B", Vec(0, 0, 0))

    simulation.registerObservable_RadialDistribution(10, rdf_callback, np.arange(0, 5, .01), "A", "B", 1. / (box_size[0] * box_size[1] * box_size[2]))
    simulation.run(T, 0.02)

    print("n_calls=%s" % n_calls)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(centers, rdf / n_calls)

    np.savetxt("bins_", centers)
    np.savetxt("rdf_", rdf/n_calls)
    plt.show()
