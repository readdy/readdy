import networkx as nx
import matplotlib.pyplot as plt


def plot_networkx_topology_graph(topology):
    G = nx.Graph()
    top_graph = topology.get_graph()
    labels = {}
    types = {}
    n_types = 0
    colors = []
    for v in top_graph.get_vertices():
        G.add_node(v.particle_index)
        labels[v.particle_index] = v.label
        if not v.particle_type() in types:
            types[v.particle_type()] = n_types
            n_types += 1
        colors.append(types[v.particle_type()])
    for v in top_graph.get_vertices():
        for vv in v:
            G.add_edge(v.particle_index, vv.get().particle_index)

    pos = nx.spring_layout(G)  # positions for all nodes

    nx.draw_networkx_nodes(G, pos, node_size=700, node_color=colors, cmap=plt.cm.summer)
    nx.draw_networkx_edges(G, pos, width=3)
    nx.draw_networkx_labels(G, pos, font_size=20, labels=labels, font_family='sans-serif')
    plt.show()
