import numpy as np
import readdy
import custom_integrator_example as cie

system = readdy.ReactionDiffusionSystem([10, 10, 10])
system.add_species("A", 0.001)
custom_integrator = cie.EBDIntegrator(.1)
sim = system.simulation(kernel="CPU")
sim.integrator = custom_integrator

sim.add_particles("A", np.random.random((3000, 3)))
sim.run(10000, .1)
