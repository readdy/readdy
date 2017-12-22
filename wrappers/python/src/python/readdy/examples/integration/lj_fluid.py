import os
import numpy as np
import matplotlib.pyplot as plt
import readdy

ut = readdy.units


def average_across_first_axis(values):
    n_values = len(values)
    mean = np.sum(values, axis=0) / n_values  # shape = n_bins
    difference = values - mean  # broadcasting starts with last axis
    std_dev = np.sqrt(np.sum(difference * difference, axis=0) / n_values)
    std_err = np.sqrt(np.sum(difference * difference, axis=0) / n_values ** 2)
    return mean, std_dev, std_err


def convert_to_numpy(vlen_array):
    nt = len(vlen_array)
    nn = len(vlen_array[0])
    result = []
    for t in range(nt):
        tmp = []
        for i in range(nn):
            tmp.append([vlen_array[t][i]['x'], vlen_array[t][i]['y'], vlen_array[t][i]['z']])
        result.append(tmp)
    result = np.array(result)
    return result


def lj_system(temperature=1.):
    system = readdy.ReactionDiffusionSystem(
        box_size=[edge_length, edge_length, edge_length],
        unit_system=None
    )
    system.kbt = temperature
    system.add_species("A", diffusion_constant=temperature)
    system.potentials.add_lennard_jones("A", "A", m=12, n=6, epsilon=1., sigma=1., cutoff=4., shift=True)
    return system


if __name__ == '__main__':
    density = 0.3
    n_particles = 1001
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


    system = lj_system(temperature=3.)

    # equilibration
    sim = system.simulation(kernel="SingleCPU")
    sim.add_particles("A", positions)

    sim.observe.particle_positions(200, callback=pos_callback, save=None)
    sim.observe.energy(50, callback=lambda x: print(x), save=None)

    sim.record_trajectory(stride=1)
    sim.output_file = "lj_eq.h5"
    if os.path.exists(sim.output_file):
        os.remove(sim.output_file)

    sim.run(n_steps=1000, timestep=2e-4)

    traj = readdy.Trajectory(sim.output_file)
    traj.convert_to_xyz(particle_radii={"A": 0.5})

    # measure
    sim = system.simulation(kernel="SingleCPU")
    sim.add_particles("A", positions)
    sim.observe.energy(50)
    sim.observe.forces(50)
    sim.observe.particle_positions(50)
    sim.observe.rdf(
        50, bin_borders=np.linspace(0.5, 4., 40),
        types_count_from="A", types_count_to="A", particle_to_density=density)

    sim.output_file = "lj_measure.h5"
    if os.path.exists(sim.output_file):
        os.remove(sim.output_file)

    sim.run(n_steps=1000, timestep=1e-4)

    traj = readdy.Trajectory(sim.output_file)
    _, energy = traj.read_observable_energy()
    _, forces = traj.read_observable_forces()
    _, positions = traj.read_observable_particle_positions()
    _, bin_centers, rdf = traj.read_observable_rdf()

    positions = convert_to_numpy(positions)
    forces = convert_to_numpy(forces)

    virial = np.sum(positions * forces, axis=2)  # sum along coordinates
    # virial = np.abs(virial)
    print("virial", virial)
    virial = np.sum(virial, axis=1) / n_particles  # average along particles
    print("virial", virial)
    mean_virial, dev_virial, err_virial = average_across_first_axis(virial)  # time average
    print("mean virial {}, err virial {}".format(mean_virial, err_virial))
    print("kT N / V", 3. * n_particles / edge_length ** 3)
    mean_pressure = 3. * n_particles / edge_length ** 3 + 1. / edge_length ** 3 * np.log(1. + 1. / 3. * mean_virial)
    err_pressure = err_virial / 3. / edge_length ** 3
    print("mean pressure {}, err pressure {}".format(mean_pressure, err_pressure))

    mean_energy, dev_energy, err_energy = average_across_first_axis(energy)  # time average
    mean_energy_per_n = mean_energy / n_particles
    dev_energy_per_n = dev_energy / n_particles
    err_energy_per_n = err_energy / n_particles
    print("mean energy per particle {}, err energy per particle {}".format(mean_energy_per_n, err_energy_per_n))

    mean_rdf, _, err_rdf = average_across_first_axis(rdf)  # time average
    plt.plot(bin_centers, mean_rdf)
    plt.show()

    mean_position = np.sum(positions, axis=1) / n_particles
    mean_position, dev_position, err_position = average_across_first_axis(mean_position)
    print("mean_position {}, dev_position {}, err_position {}".format(mean_position, dev_position, err_position))
