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

        recorded_topologies = []
        recorded_particles = []

        sim.make_checkpoints(10)
        sim.record_trajectory()
        sim.observe.topologies(1, callback=lambda x: recorded_topologies.append(x))
        sim.observe.particles(1, callback=lambda x: recorded_particles.append(x))
        sim.output_file = os.path.join(self.dir, 'traj.h5')
        sim.show_progress = False
        sim.run(500, 1e-3, show_system=False)

        checkpoint_step = sim.list_checkpoints(sim.output_file)[3][1]['step']

        system = self._set_up_system()
        sim = system.simulation()
        sim.load_particles_from_checkpoint(os.path.join(self.dir, 'traj.h5'), 3)

        tops = sim.current_topologies
        assert len(tops) == 2
        assert tops[0].type == "TT1"
        assert tops[0].graph.has_edge(0, 1)
        assert tops[0].graph.has_edge(1, 2)
        assert tops[0].graph.has_edge(2, 3)
        assert tops[0].graph.has_edge(3, 0)
        assert not tops[0].graph.has_edge(0, 2)
        assert not tops[0].graph.has_edge(1, 3)
        assert tops[1].type == "TT2"
        assert tops[1].graph.has_edge(0, 1)
        assert tops[1].graph.has_edge(1, 2)
        assert tops[1].graph.has_edge(2, 3)

        p_types, p_ids, p_positions = recorded_particles[checkpoint_step]
        topologies = recorded_topologies[checkpoint_step]

        assert len(topologies) == len(tops)
        for ix, topology_record in enumerate(topologies):
            for edge in topology_record.edges:
                assert tops[ix].graph.has_edge(*edge)
            for particle in topology_record.particles:
                top_particle = next(p for p in tops[ix].particles if p.id == p_ids[particle])
                assert top_particle.id == p_ids[particle]
                assert top_particle.type == p_types[particle]
                assert top_particle.pos == p_positions[particle]


        sim.show_progress = False
        sim.run(500, 1e-3, show_system=False)
