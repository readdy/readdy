from readdy._internal.simulation import KernelProvider, Simulation, Vec

from readdy.util import platform_utils

if __name__ == '__main__':
    kernel_provider = KernelProvider.get()
    kernel_provider.load_from_dir(platform_utils.get_readdy_plugin_dir())
    simulation = Simulation()
    simulation.setKernel("SingleCPU")

    box_size = Vec(20, 20, 20)
    # membrane particle
    simulation.registerParticleType("M", 0, 1)
    # MinD-ADP (without phosphor)
    simulation.registerParticleType("D", 2.5, .5)
    # MinD-ATP (with phosphor)
    simulation.registerParticleType("D_P", 2.5, .5)
    # MinE
    simulation.registerParticleType("E", 2.5, .5)
    # MinD-ATP bound
    simulation.registerParticleType("D_PB", .01, .5)
    # MinDE
    simulation.registerParticleType("DE", .01, .5)

    simulation.registerConversionReaction("D", "D_P", .5)