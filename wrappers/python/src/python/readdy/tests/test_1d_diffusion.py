import unittest
import numpy as np
import readdy


class Test1DDiffusion(unittest.TestCase):
    def test_1d_diffusion(self):
        sys = readdy.ReactionDiffusionSystem([10, 10, 10])
        sys.add_species("A", [0., 0., 1.])
        sim = sys.simulation()
        sim.add_particle("A", [0., 0., 0.])

        traj = []
        sim.observe.particle_positions(1, callback=lambda x: traj.append(x))

        sim.run(100, 1e-3, show_summary=False)

        for x in traj:
            np.testing.assert_equal(x[0][0], 0)
            np.testing.assert_equal(x[0][1], 0)
