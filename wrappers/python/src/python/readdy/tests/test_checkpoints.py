import os
import shutil
import tempfile

import h5py

import readdy
from readdy.util.testing_utils import ReaDDyTestCase

import numpy as np

class TestCheckpoints(ReaDDyTestCase):

    @classmethod
    def setUpClass(cls):
        cls.dir = tempfile.mkdtemp("test-checkpoints")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.dir, ignore_errors=True)

    def _set_up_system(self):
        system = readdy.ReactionDiffusionSystem(box_size=[10, 10, 10])
        system.add_species("A", 1.)
        system.add_species("B", 1.)
        system.topologies.add_type("TT1")
        system.topologies.add_type("TT2")
        system.add_topology_species("T1", 1.)
        system.add_topology_species("T2", 1.)
        system.topologies.configure_harmonic_bond("T1", "T1")
        system.topologies.configure_harmonic_bond("T1", "T2")
        system.topologies.configure_harmonic_bond("T2", "T2")
        return system

    def test_continue_simulation_no_free_particles(self):
        pass

    def test_continue_simulation_no_topologies(self):
        pass

    def test_continue_simulation_full(self):
        system = self._set_up_system()

        print(system.kbt)

        sim = system.simulation()

        t1 = sim.add_topology("TT1", ["T1", "T2", "T1", "T2"], np.random.normal(0, 1, size=(4, 3)))
        t1.graph.add_edge(0, 1)
        t1.graph.add_edge(1, 2)
        t1.graph.add_edge(2, 3)
        t1.graph.add_edge(3, 0)
        t2 = sim.add_topology("TT2", ["T2", "T1", "T2", "T1"], np.random.normal(0, 1, size=(4, 3)))
        t2.graph.add_edge(0, 1)
        t2.graph.add_edge(1, 2)
        t2.graph.add_edge(2, 3)

        sim.add_particles("A", np.random.normal(0, 1, size=(20, 3)))
        sim.add_particles("B", np.random.normal(0, 1, size=(50, 3)))

        sim.make_checkpoints(100)
        sim.record_trajectory()
        sim.observe.topologies(1)
        sim.output_file = os.path.join(self.dir, 'traj.h5')
        sim.run(500, 1e-3, show_system=False)

        def visitor_func(name, node):
            if isinstance(node, h5py.Dataset):
                print(f"\tData set {name}")
            else:
                # node is a group
                print(f"\tGroup {name}")

        f = h5py.File(sim.output_file,'r')
        f.visititems(visitor_func)

        print(sim.list_checkpoints(sim.output_file))
