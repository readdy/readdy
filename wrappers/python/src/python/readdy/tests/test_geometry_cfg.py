import unittest
import readdy
import numpy as np

class TestGeometryInterface(unittest.TestCase):

    def test_box_potential(self):
        system = readdy.ReactionDiffusionSystem([10, 10, 10])
        system.add_species("A", diffusion_constant=1.)
        system.potentials.add_harmonic_geometry("A", 100, system.geometry.create_box([-1, -1, -1], [3, 3, 3]), True)
        system.potentials.add_harmonic_geometry("A", 100, system.geometry.create_sphere([0, 0, 0], 1.), True)
        system.potentials.add_harmonic_geometry("A", 100, system.geometry.create_capsule([0., 0., 0.], [1., 0., 0.], 1., 2.), True)
        sim = system.simulation()
        sim.add_particles("A", np.random.uniform(-5, 5, size=(100, 3)))
        positions = []
        sim.observe.particle_positions(stride=100, callback=lambda posns: positions.append(posns))
        sim.run(10000, timestep=1e-4)

    def test_box_compartment(self):
        system = readdy.ReactionDiffusionSystem([10, 10, 10])
        system.add_species("A", diffusion_constant=1.)
        system.add_species("Aint", diffusion_constant=1.)
        system.compartments.add_geometry({"A": "Aint"}, "tmp", system.geometry.create_box([-1, -1, -1], [3, 3, 3]), True)
        system.compartments.add_geometry({"A": "Aint"}, "tmp", system.geometry.create_sphere([0, 0, 0], 1.), True)
        system.compartments.add_geometry({"A": "Aint"}, "tmp", system.geometry.create_capsule([0., 0., 0.], [1., 0., 0.], 1., 2.), True)
        sim = system.simulation()
        sim.add_particles("A", np.random.uniform(-5, 5, size=(100, 3)))
        sim.run(1000, timestep=1e-3)
