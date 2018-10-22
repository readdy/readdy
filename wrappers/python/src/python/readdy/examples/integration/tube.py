import os
import numpy as np
import readdy

rds = readdy.ReactionDiffusionSystem(box_size=[10., 10., 10.], unit_system=None)
rds.add_species("A")
rds.potentials.add_cylinder("A", 100, (0, 0, 0), (1, 0, 0), 1, inclusion=False)

simulation = rds.simulation(kernel="SingleCPU")

simulation.output_file = "out.h5"
simulation.add_particles("A", np.random.normal(size=(5000, 3)))
simulation.record_trajectory(stride=100)

if os.path.exists(simulation.output_file):
    os.remove(simulation.output_file)
simulation.run(n_steps=30000, timestep=1e-3)

trajectory = readdy.Trajectory('out.h5')
trajectory.convert_to_xyz(particle_radii={'A': 0.1})
