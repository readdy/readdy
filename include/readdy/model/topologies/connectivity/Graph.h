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
#include <unordered_map>
#include <readdy/common/macros.h>
#include "Vertex.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(top)
NAMESPACE_BEGIN(graph)

class Graph {
public:
    Graph() = default;

    virtual ~Graph() = default;

    Graph(const Graph &) = delete;

    Graph &operator=(const Graph &) = delete;

    Graph(Graph&&) = default;

    Graph &operator=(Graph&&) = default;

    const std::vector<Vertex> &vertices() const;

    const Vertex* const namedVertex(const std::string& name) const;

    const Vertex* const vertexForParticleIndex(std::size_t particleIndex) const;

    void addVertex(const Vertex &);

    void addVertex(Vertex &&);

    void setVertexLabel(std::size_t vertex, const std::string& label);

    void addEdge(std::size_t v1, std::size_t v2);

    void addEdge(std::size_t v1, const std::string &v2);

    void addEdge(const std::string &v1, std::size_t v2);

    void addEdge(const std::string &v1, const std::string &v2);

    void addEdgeBetweenParticles(std::size_t particleIndex1, std::size_t particleIndex2);

    void removeEdge(std::size_t v1, std::size_t v2);

    void removeEdge(std::size_t v1, const std::string &v2);

    void removeEdge(const std::string &v1, std::size_t v2);

    void removeEdge(const std::string &v1, const std::string &v2);

    void removeVertex(std::size_t index);

    void removeVertex(Vertex *vertex);

    void removeVertex(const std::string &name);

    void removeParticle(std::size_t particleIndex);

private:
    std::vector<Vertex> vertices_;
    std::unordered_map<std::string, Vertex *> namedVertices;

    void removeNeighborsEdges(Vertex *vertex);

    void removeEdge(Vertex *v1, Vertex *v2);

    auto vertexItForParticleIndex(std::size_t particleIndex) -> decltype(vertices_.begin());
};

NAMESPACE_END(graph)
NAMESPACE_END(top)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
