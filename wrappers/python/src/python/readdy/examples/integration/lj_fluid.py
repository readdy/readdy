import os
import numpy as np
import sys
sys.path.append("/home/chris/workspace/readdy/wrappers/python/src/python")

import readdy

ut = readdy.units


def lj_system(temperature=1):
    system = readdy.ReactionDiffusionSystem(
        box_size=[edge_length, edge_length, edge_length],
        unit_system=None
    )
    system.kbt = temperature
    system.add_species("A", diffusion_constant=1.)
    system.potentials.add_lennard_jones("A", "A", m=12, n=6, epsilon=1., sigma=1., cutoff=None, shift=True)
    return system


density = 0.3
n_particles = 1000
edge_length = (n_particles / density) ** (1. / 3.)

n_per_axis = int(n_particles ** (1. / 3.))
pos_x = np.linspace(-edge_length / 2., edge_length / 2. - 1., n_per_axis)
print(n_per_axis)
positions = []
for x in pos_x:
    for y in pos_x:
        for z in pos_x:
            positions.append([x, y, z])
positions = np.array(positions)


def pos_callback(x):
    global positions
    n = len(x)
    positions = np.zeros((n, 3))
    for i in range(n):
        positions[i][0] = x[i][0]
        positions[i][1] = x[i][1]
        positions[i][2] = x[i][2]
    print("saved positions")


system = lj_system(3.)
sim = system.simulation(kernel="SingleCPU")
sim.add_particles("A", positions)

sim.observe.particle_positions(5000, callback=pos_callback, save=None)
sim.observe.energy(500, callback=lambda x: print(x), save=None)

sim.record_trajectory(stride=1)
sim.output_file = "lj_eq.h5"
if os.path.exists(sim.output_file):
    os.remove(sim.output_file)

sim.run(n_steps=10000, timestep=2e-5)

traj = readdy.Trajectory(sim.output_file)
traj.convert_to_xyz(particle_radii={"A": 0.5})
