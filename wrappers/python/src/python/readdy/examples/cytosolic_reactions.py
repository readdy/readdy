# coding=utf-8

# Copyright © 2018 Computational Molecular Biology Group,
#                  Freie Universität Berlin (GER)
#
# Redistribution and use in source and binary forms, with or
# without modification, are permitted provided that the
# following conditions are met:
#  1. Redistributions of source code must retain the above
#     copyright notice, this list of conditions and the
#     following disclaimer.
#  2. Redistributions in binary form must reproduce the above
#     copyright notice, this list of conditions and the following
#     disclaimer in the documentation and/or other materials
#     provided with the distribution.
#  3. Neither the name of the copyright holder nor the names of
#     its contributors may be used to endorse or promote products
#     derived from this software without specific
#     prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
# CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
# INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
# STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
# ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import os
import socket
import argparse
import pickle
import numpy as np
import time
import multiprocessing

import readdy

"""
Usage: `python cytosolic_reactions.py $version_string`

Perform simulations for system A + B <--> C and repulsion between all particles.
This is done for different number of particles (at fixed density) and different compute-kernels.
The benchmark value for the different compute-kernels is the computation-time per particle and per integration step
at the highest number of particles (where this benchmark value is typically best, due to load balancing).

Output is a dictionary that contains relevant performance data for each simulation as well as the version string and
the name of the host machine. The output is pickled to compare different versions later.
"""

parser = argparse.ArgumentParser(description='Run performance scenario A + B <--> C with repulsion')
parser.add_argument('--long', action='store_true',
                    help='Perform full equilibration of the system. Computes for several hours')
parser.add_argument('--debug', action='store_true',
                    help='Very short run for debugging output file and stuff. Computes for few minutes')


def traverse_performance_tree(tree):
    if len(tree) == 0:
        # leaf
        return {"time": tree.time(), "count": tree.count()}
    else:
        result = {"children": {}, "time": tree.time(), "count": tree.count()}
        for name in tree.keys():
            result["children"][name] = traverse_performance_tree(tree[name])
        return result


def perform(kernel="SingleCPU", n_particles_a=2357, force_constant=10., file_suffix="", full_simulation=False,
            debug_run=False, n_threads=-1):
    print("kernel {}, n_particles_a {}, force_constant {}, threads {}"
          .format(kernel, n_particles_a, force_constant, n_threads))
    n_particles_b = n_particles_a
    n_particles_c = 0
    desired_a_density = 2357. / 1e6  # number of particles per nanometer**3
    edge_length = (n_particles_a / desired_a_density) ** (1. / 3.)
    print("Use edge length {}".format(edge_length))
    system = readdy.ReactionDiffusionSystem(
        box_size=[edge_length, edge_length, edge_length],
        temperature=293.,
        unit_system={"length_unit": "nanometer", "time_unit": "nanosecond", "energy_unit": "kilojoule / mol"}
    )

    particle_radii = {"A": 1.5, "B": 3., "C": 3.12}  # in nanometers
    ut = readdy.units

    system.add_species("A", diffusion_constant=143.1 * ut.micrometer ** 2 / ut.second)
    system.add_species("B", diffusion_constant=71.6 * ut.micrometer ** 2 / ut.second)
    system.add_species("C", diffusion_constant=68.82 * ut.micrometer ** 2 / ut.second)

    reaction_radius = particle_radii["A"] + particle_radii["B"]
    system.reactions.add_fusion("fusion", "A", "B", "C",
                                rate=1e6 * 1. / ut.second,
                                educt_distance=reaction_radius * ut.nanometer)
    system.reactions.add_fission("fission", "C", "A", "B",
                                 rate=5e4 * 1. / ut.second,
                                 product_distance=reaction_radius * ut.nanometer)

    if force_constant > 0.:
        fc = force_constant * ut.kilojoule / ut.mol / (ut.nanometer ** 2)
        system.potentials.add_harmonic_repulsion("A", "A", fc,
                                                 interaction_distance=particle_radii["A"] + particle_radii["A"])
        system.potentials.add_harmonic_repulsion("B", "B", fc,
                                                 interaction_distance=particle_radii["B"] + particle_radii["B"])
        system.potentials.add_harmonic_repulsion("C", "C", fc,
                                                 interaction_distance=particle_radii["C"] + particle_radii["C"])
        system.potentials.add_harmonic_repulsion("A", "B", fc,
                                                 interaction_distance=particle_radii["A"] + particle_radii["B"])
        system.potentials.add_harmonic_repulsion("B", "C", fc,
                                                 interaction_distance=particle_radii["B"] + particle_radii["C"])
        system.potentials.add_harmonic_repulsion("C", "A", fc,
                                                 interaction_distance=particle_radii["C"] + particle_radii["A"])

    simulation = system.simulation(kernel=kernel)
    simulation.output_file = "cytosolic_reactions_" + kernel + "_n_a_" + str(n_particles_a) \
                             + "_force_" + str(force_constant) + "_" + file_suffix + ".h5"
    simulation.reaction_handler = "UncontrolledApproximation"

    edge_length = system.box_size[0]
    if full_simulation:
        initial_positions_a = np.random.random(size=(n_particles_a, 3)) * edge_length - .5 * edge_length
        initial_positions_b = np.random.random(size=(n_particles_b, 3)) * edge_length - .5 * edge_length
        simulation.add_particles("A", initial_positions_a)
        simulation.add_particles("B", initial_positions_b)
    else:  # start in roughly equilibrated state
        n_particles_c = n_particles_a * 2 // 3
        n_particles_a_actual = n_particles_a - n_particles_c
        n_particles_b = n_particles_a_actual
        initial_positions_a = np.random.random(size=(n_particles_a_actual, 3)) * edge_length - .5 * edge_length
        initial_positions_b = np.random.random(size=(n_particles_b, 3)) * edge_length - .5 * edge_length
        initial_positions_c = np.random.random(size=(n_particles_c, 3)) * edge_length - .5 * edge_length
        simulation.add_particles("A", initial_positions_a)
        simulation.add_particles("B", initial_positions_b)
        simulation.add_particles("C", initial_positions_c)

    simulation.observe.number_of_particles(stride=1, types=["A", "B", "C"])

    if os.path.exists(simulation.output_file):
        os.remove(simulation.output_file)

    dt = 1e-1  # in nanoseconds
    if full_simulation:
        n_steps = int(10000. / dt)  # simulate to 10 microseconds
    else:
        n_steps = 3000
    if debug_run:
        n_steps = 200
    print("Performing n_steps {} ..".format(n_steps))

    if kernel != 'SingleCPU':
        simulation.kernel_configuration.n_threads = n_threads

    simulation.run(n_steps=n_steps, timestep=dt * ut.nanosecond)

    perf = simulation._simulation.performance_root()
    print(perf)
    print("Simulated for {} seconds".format(perf.time()))
    performance_tree = traverse_performance_tree(perf)

    traj = readdy.Trajectory(simulation.output_file)
    times, counts = traj.read_observable_number_of_particles()
    counts = {"A": counts[:, 0], "B": counts[:, 1], "C": counts[:, 2]}
    times = times * dt
    n_frames = len(counts["A"])
    average_n_particles = (np.sum(counts["A"]) + np.sum(counts["B"]) + np.sum(counts["C"])) / n_frames
    print("Time averaged total number of particles {}".format(average_n_particles))

    t_pp_ps = perf.time() / average_n_particles / n_steps
    print("Computation time per particle per integration step is {} seconds".format(t_pp_ps))

    arguments = {
        "kernel": kernel,
        "n_particles_a": n_particles_a,
        "force_constant": force_constant,
        "file_suffix": file_suffix
    }
    result = {
        "time/particle/step": t_pp_ps,
        "average_n_particles": average_n_particles,
        "n_particles": counts,
        "times": times,
        "computation_time": perf.time(),
        "performance_tree": performance_tree,
        "volume": edge_length ** 3,
        "unit_system": {"length_unit": "nanometer", "time_unit": "nanosecond", "energy_unit": "kilojoule / mol"},
        "n_steps": n_steps,
        "n_threads": n_threads
    }
    data = dict()
    data.update(arguments)
    data.update(result)

    os.unlink(simulation.output_file)

    del simulation
    del system

    import gc
    gc.collect()
    return data


if __name__ == '__main__':
    print("run cytosolic_reactions ..")
    args = parser.parse_args()
    full_simulation = args.long
    debug_run = args.debug
    version = readdy.__version__
    host = socket.gethostname()
    file_suffix = version + "_" + host
    file_name = "cytosolic_reactions_" + file_suffix + ".pickle"
    if os.path.exists(file_name):
        raise RuntimeError("An output file with the desired name {} already exists. ".format(file_name))

    t1 = time.perf_counter()
    n_cores = multiprocessing.cpu_count()
    perf_results = []
    for kernel in ["CPU", "SingleCPU"]:
            for n_particles_a in [200, 400, 700, 1000, 2357, 5000, 10000, 15000, 20000, 30000]:
                for force_constant in [10.]:
                    if kernel != "SingleCPU":
                        for n_threads in [1,2,3,4,5,6,7,8,10,16,24,32]:
                            data = perform(
                                kernel=kernel,
                                n_particles_a=n_particles_a,
                                force_constant=force_constant,
                                file_suffix=version,
                                full_simulation=full_simulation,
                                debug_run=debug_run,
                                n_threads=n_threads)
                            perf_results.append(data)
                            time.sleep(.5)
                    else:
                        data = perform(
                            kernel=kernel,
                            n_particles_a=n_particles_a,
                            force_constant=force_constant,
                            file_suffix=version,
                            full_simulation=full_simulation,
                            debug_run=debug_run)
                        perf_results.append(data)

    t2 = time.perf_counter()
    print("The measuring procedure took {} seconds".format(t2 - t1))

    file_data = {"results": perf_results, "version": version, "host": host}
    with open(file_name, "xb") as f:
        pickle.dump(file_data, f, pickle.HIGHEST_PROTOCOL)