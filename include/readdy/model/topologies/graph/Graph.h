/********************************************************************
 * Copyright © 2016 Computational Molecular Biology Group,          * 
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * This file is part of ReaDDy.                                     *
 *                                                                  *
 * ReaDDy is free software: you can redistribute it and/or modify   *
 * it under the terms of the GNU Lesser General Public License as   *
 * published by the Free Software Foundation, either version 3 of   *
 * the License, or (at your option) any later version.              *
 *                                                                  *
 * This program is distributed in the hope that it will be useful,  *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of   *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    *
 * GNU Lesser General Public License for more details.              *
 *                                                                  *
 * You should have received a copy of the GNU Lesser General        *
 * Public License along with this program. If not, see              *
 * <http://www.gnu.org/licenses/>.                                  *
 ********************************************************************/


/**
 * << detailed description >>
 *
 * @file Graph.h
 * @brief << brief description >>
 * @author clonker
 * @date 16.03.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <functional>
#include <stdexcept>
#include <list>
#include <unordered_map>
#include <readdy/common/macros.h>
#include "Vertex.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(top)
NAMESPACE_BEGIN(graph)

class Graph {
public:

    using vertex_list = std::list<Vertex>;
    using vertex_ref = vertex_list::iterator;
    using vertex_cref = vertex_list::const_iterator;

    using edge = std::tuple<vertex_ref, vertex_ref>;
    using cedge = std::tuple<vertex_cref, vertex_cref>;

    using path_len_2 = std::tuple<vertex_ref, vertex_ref, vertex_ref>;
    using cpath_len_2 = std::tuple<vertex_cref, vertex_cref, vertex_cref>;

    using path_len_3 = std::tuple<vertex_ref, vertex_ref, vertex_ref, vertex_ref>;
    using cpath_len_3 = std::tuple<vertex_cref, vertex_cref, vertex_cref, vertex_cref>;

    using edge_callback = std::function<void(const edge &)>;
    using path_len_2_callback = std::function<void(const path_len_2 &)>;
    using path_len_3_callback = std::function<void(const path_len_3 &)>;

    Graph() = default;

    explicit Graph(vertex_list vertexList);

    virtual ~Graph() = default;

    Graph(const Graph &) = delete;

    Graph &operator=(const Graph &) = delete;

    Graph(Graph &&) = default;

    Graph &operator=(Graph &&) = default;

    const vertex_list &vertices() const;

    vertex_list &vertices();

    vertex_ref firstVertex();

    vertex_ref lastVertex();

    bool containsEdge(const cedge& edge) const;

    bool containsEdge(vertex_cref v1, vertex_cref v2) const;

    const Vertex &vertexForParticleIndex(std::size_t particleIndex) const;

    void addVertex(std::size_t particleIndex, particle_type_type particleType);

    void addEdge(vertex_ref v1, vertex_ref v2);

    void addEdge(const edge& edge);

    void addEdgeBetweenParticles(std::size_t particleIndex1, std::size_t particleIndex2);

    void removeEdge(vertex_ref v1, vertex_ref v2);

    void removeEdge(const edge& edge);

    void removeVertex(vertex_ref vertex);

    void removeParticle(std::size_t particleIndex);

    bool isConnected();
    
    std::vector<std::tuple<vertex_ref, vertex_ref>> edges();

    void findEdges(const edge_callback &edgeCallback);
    
    void findNTuples(const edge_callback &tuple_callback,
                     const path_len_2_callback &triple_callback,
                     const path_len_3_callback &quadruple_callback);

    std::tuple<std::vector<edge>, std::vector<path_len_2>, std::vector<path_len_3>>
    findNTuples();

    /**
     * Returns the connected components, invalidates this graph
     * @return connected components
     */
    std::vector<Graph> connectedComponentsDestructive();

private:
    vertex_list _vertices {};

    void removeNeighborsEdges(vertex_ref vertex);

    auto vertexItForParticleIndex(std::size_t particleIndex) -> decltype(_vertices.begin());
};

NAMESPACE_END(graph)
NAMESPACE_END(top)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
