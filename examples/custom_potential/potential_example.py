import numpy as np
import readdy
import custom_potential_example as cb
system = readdy.ReactionDiffusionSystem([10, 10, 10])
system.add_species("A", 0.001)
system.potentials.add_custom_external("A", cb.TestPotential)
system.potentials.add_custom_pair("A", "A", cb.TestPairPotential, .1)
sim = system.simulation(kernel="CPU")

sim.add_particles("A", np.random.random((3000, 3)))
sim.run(10000, .1)
