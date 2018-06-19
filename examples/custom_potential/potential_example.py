import numpy as np
import readdy
import cudakernel_binding as cb
system = readdy.ReactionDiffusionSystem([10, 10, 10])
system.add_species("A", 0.001)
p = cb.TestPotential(0)
system._context.potentials.add_external_order1(p)
sim = system.simulation(kernel="CPU")

sim.add_particles("A", np.random.random((30000, 3)))
sim.run(10000, .1)
