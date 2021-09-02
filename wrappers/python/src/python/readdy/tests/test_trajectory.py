import shutil
import tempfile
import unittest
import readdy
import numpy as np


class TestTrajectory(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.dir = tempfile.mkdtemp("test-observables-io")
        sys = readdy.ReactionDiffusionSystem([10, 10, 10])
        sys.add_species("A", 1.)
        sim = sys.simulation()
        sim.record_trajectory()
        sim.add_particles("A", np.random.uniform(-4, 4, size=(100, 3)))
        sim.output_file = cls.dir + "/out.h5"
        sim.run(123, 1e-4)
        cls.traj = readdy.Trajectory(sim.output_file)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.dir, ignore_errors=True)

    def test_trajectory_length(self):
        self.assertEqual(len(self.traj), 123)
