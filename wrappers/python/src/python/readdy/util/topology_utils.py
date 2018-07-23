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

def plot_networkx_topology_graph(topology):
    import matplotlib.pyplot as plt
    import networkx as nx

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


def plot_networkx_graph(G):
    import matplotlib.pyplot as plt
    import networkx as nx

    pos = nx.spring_layout(G)  # positions for all nodes
    labels = {}
    for node in G.nodes():
        labels[node] = G.node[node]["label"]
        if not labels[node]:
            labels[node] = node
    nx.draw_networkx_nodes(G, pos, node_size=700, cmap=plt.cm.summer)
    nx.draw_networkx_edges(G, pos, width=3)
    nx.draw_networkx_labels(G, pos, font_size=20, labels=labels, font_family='sans-serif')
    plt.show()


def plot_gexf_string(string):
    import networkx as nx
    from io import StringIO

    strio = StringIO(u"%s" % string)
    graph = nx.read_gexf(strio, relabel=False)
    plot_networkx_graph(graph)
