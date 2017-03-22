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
    using vertices_t = std::list<Vertex>;

    Graph() = default;

    virtual ~Graph() = default;

    Graph(const Graph &) = delete;

    Graph &operator=(const Graph &) = delete;

    Graph(Graph &&) = default;

    Graph &operator=(Graph &&) = default;

    const vertices_t &vertices() const;

    vertices_t &vertices();

    vertices_t::iterator firstVertex();

    vertices_t::iterator lastVertex();

    const Vertex &namedVertex(const std::string &name) const;

    Vertex &namedVertex(const std::string &name);

    const Vertex &vertexForParticleIndex(std::size_t particleIndex) const;

    void addVertex(std::size_t particleIndex, particle_type_type particleType, const std::string &label = "");

    void setVertexLabel(vertices_t::iterator vertex, const std::string &label);

    void addEdge(vertices_t::iterator v1, vertices_t::iterator v2);

    void addEdge(const std::string &v1, const std::string &v2);

    void addEdgeBetweenParticles(std::size_t particleIndex1, std::size_t particleIndex2);

    void removeEdge(vertices_t::iterator v1, vertices_t::iterator v2);

    void removeEdge(const std::string &v1, const std::string &v2);

    void removeVertex(vertices_t::iterator vertex);

    void removeVertex(const std::string &name);

    void removeParticle(std::size_t particleIndex);

    bool isConnected();

private:
    vertices_t vertices_;
    std::unordered_map<std::string, vertices_t::iterator> namedVertices{};

    void removeNeighborsEdges(vertices_t::iterator vertex);

    auto vertexItForParticleIndex(std::size_t particleIndex) -> decltype(vertices_.begin());
};

NAMESPACE_END(graph)
NAMESPACE_END(top)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
