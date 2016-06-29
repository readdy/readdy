from readdy._internal.simulation import KernelProvider, Simulation, Vec

from readdy.util import platform_utils


class MinEMinDSimulation(object):

    def com_callback(self, centerOfMass):
        print("center of mass = %s" % centerOfMass)

    def execute(self):
        kernel_provider = KernelProvider.get()
        kernel_provider.load_from_dir(platform_utils.get_readdy_plugin_dir())
        simulation = Simulation()
        simulation.setKernel("SingleCPU")

        ###################################
        #
        # set up simulation box
        #
        ###################################

        box_size = Vec(20, 20, 20)
        simulation.box_size = box_size
        simulation.kbt = 1  # todo: resonable?
        simulation.periodic_boundary = [False, False, False]

        ###################################
        #
        # register particle types
        #
        ###################################

        simulation.registerParticleType("M", 0, 1)  # membrane particle
        simulation.registerParticleType("D", 2.5, .5)  # MinD-ADP (without phosphor)
        simulation.registerParticleType("D_P", 2.5, .5)  # MinD-ATP (with phosphor)
        simulation.registerParticleType("E", 2.5, .5)  # MinE
        simulation.registerParticleType("D_PB", .01, .5)  # MinD-ATP bound
        simulation.registerParticleType("DE", .01, .5)  # MinDE

        ###################################
        #
        # register reaction types
        #
        ###################################

        simulation.registerConversionReaction("Phosphorization", "D", "D_P", .5)
        simulation.registerEnzymaticReaction("Attach to membrane", "M", "D_P", "D_PB", .025, .25)  # todo: check rate, check educt distance
        simulation.registerFusionReaction("bound MinD+MinE->MinDE", "D_PB", "E", "DE", .093, .25, .5, .5)  # todo: apply erban, chapman to rate, educt distance?
        simulation.registerFissionReaction("MinDE to MinD and MinE, detach", "DE", "D", "E", .25, .7, .5, .5)  # todo: reasonable product distance?

        ###################################
        #
        # register observables
        #
        ###################################

        simulation.registerObservable_CenterOfMass(1, self.com_callback, ["D", "D_P", "D_PB"])

        ###################################
        #
        # register potentials
        #
        ###################################


if __name__ == '__main__':
    sim = MinEMinDSimulation()
    sim.execute()
