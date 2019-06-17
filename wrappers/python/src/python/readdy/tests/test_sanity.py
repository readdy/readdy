import os
import shutil
import tempfile
import unittest

import readdy


class SanityTests(unittest.TestCase):

    @classmethod
    def setUpClass(cls) -> None:
        cls.dir = tempfile.mkdtemp("test-sanity")

    @classmethod
    def tearDownClass(cls) -> None:
        shutil.rmtree(cls.dir, ignore_errors=True)

    def test_run_without_progress(self):
        # deals with issue 161
        system = readdy.ReactionDiffusionSystem([25.,25.,25.], temperature=300.*readdy.units.kelvin)
        simulation = system.simulation(kernel="CPU")
        simulation.output_file = os.path.join(self.dir, 'test_run_without_progress.h5')

        simulation.show_progress = False
        simulation.run(n_steps=10000, timestep=1e-2)
