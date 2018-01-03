import readdy

if __name__ == '__main__':
    import numpy as np

    rds = readdy.ReactionDiffusionSystem(box_size=[10., 10., 10.])
    rds.add_species("A", diffusion_constant=1.0)
    rds.add_species("B", diffusion_constant=1.0)
    rds.add_topology_species("T", diffusion_constant=5)
    rds.topologies.add_type("TT")
    rds.reactions.add_conversion("myconversion", "A", "B", 1.0)
    sim = rds.simulation(kernel="CPU")
    t = sim.add_topology("TT", ["T", "T", "T"], np.random.random((3, 3)))
    graph = t.get_graph()
    graph.add_edge(0, 1)
    graph.add_edge(1, 2)
    sim.run(5, .1)
