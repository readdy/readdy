import os
import shutil
import tempfile

import numpy as np

import readdy
from readdy.util.testing_utils import ReaDDyTestCase


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
        system.reactions.add("conv1: A -> B", rate=1.)
        system.reactions.add("conv2: B -> A", rate=1.)
        system.reactions.add("fus: A +(1.) B -> A", rate=1.)
        system.reactions.add("fiss: A -> A +(1.) B", rate=1.)
        system.topologies.add_type("TT1")
        system.topologies.add_type("TT2")
        system.topologies.add_type("Dummy")
        system.add_topology_species("Dummy", 10.)
        system.add_topology_species("T1", 1.)
        system.add_topology_species("T2", 1.)
        system.topologies.configure_harmonic_bond("T1", "T1")
        system.topologies.configure_harmonic_bond("T1", "T2")
        system.topologies.configure_harmonic_bond("T2", "T2")
        system.topologies.configure_harmonic_bond("Dummy", "Dummy")
        return system

    def test_continue_simulation_full(self):
        system = self._set_up_system()
        sim = system.simulation()

        t1_initial_pos = np.random.normal(0, 1, size=(4, 3))
        t1 = sim.add_topology("TT1", ["T1", "T2", "T1", "T2"], t1_initial_pos)
        t1.graph.add_edge(0, 1)
        t1.graph.add_edge(1, 2)
        t1.graph.add_edge(2, 3)
        t1.graph.add_edge(3, 0)
        t2_initial_pos = np.random.normal(0, 1, size=(4, 3))
        t2 = sim.add_topology("TT2", ["T2", "T1", "T2", "T1"], t2_initial_pos)
        t2.graph.add_edge(0, 1)
        t2.graph.add_edge(1, 2)
        t2.graph.add_edge(2, 3)

        a_particles_initial_pos = np.random.normal(0, 1, size=(20, 3))
        sim.add_particles("A", a_particles_initial_pos)
        b_particles_initial_pos = np.random.normal(0, 1, size=(50, 3))
        sim.add_particles("B", b_particles_initial_pos)

        recorded_topologies = []

        def topologies_callback(x):
            recorded_topologies.append(x)
            if len(sim.current_topologies) % 2 == 0:
                sim.add_topology("Dummy", "Dummy", np.random.random(size=(1, 3)))
            else:
                t = sim.add_topology("Dummy", "Dummy", np.random.random(size=(5, 3)))
                t.graph.add_edge(0, 1)
                t.graph.add_edge(1, 2)
                t.graph.add_edge(2, 3)
                t.graph.add_edge(3, 4)
                t.configure()

        sim.make_checkpoints(7)
        sim.record_trajectory()
        sim.observe.topologies(1, callback=topologies_callback)
        sim.output_file = os.path.join(self.dir, 'traj.h5')
        sim.show_progress = False
        sim.run(120, 1e-2, show_system=False)

        traj = readdy.Trajectory(sim.output_file)
        entries = traj.read()
        checkpoint = sim.list_checkpoints(sim.output_file)[17]

        system = self._set_up_system()
        sim = system.simulation()
        sim.load_particles_from_checkpoint(os.path.join(self.dir, 'traj.h5'), checkpoint['number'])

        tops = sim.current_topologies
        assert len(tops) == 2 + checkpoint['step'], f"expected {2 + checkpoint['step']} topologies, " \
                                                    f"got {len(tops)}"
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

        current_entries = entries[checkpoint['step']]
        current_particles = sim.current_particles
        topologies = recorded_topologies[checkpoint['step']]

        # check whether restored topologies are o.k.
        assert len(topologies) == len(tops)
        for ix, topology_record in enumerate(topologies):
            restored_topology = tops[ix]
            for edge in topology_record.edges:
                assert tops[ix].graph.has_edge(*edge)
            for pix, particle_ix in enumerate(topology_record.particles):
                particle = current_entries[particle_ix]
                restored_particle = restored_topology.particles[pix]
                assert np.array_equal(restored_particle.pos, np.array(particle.position))
                assert restored_particle.type == particle.type

        # check whether restored particles are o.k.
        for entry in current_entries:
            # see if entry is available in current particles
            ix = 0
            for ix, particle in enumerate(current_particles):
                if particle.type == entry.type and np.array_equal(particle.pos, entry.position):
                    break
            assert ix < len(current_particles), f"entry {entry} was not found in particles!"
            current_particles.pop(ix)
        assert len(current_particles) == 0

        sim.show_progress = False
        sim.run(500, 1e-3, show_system=False)
