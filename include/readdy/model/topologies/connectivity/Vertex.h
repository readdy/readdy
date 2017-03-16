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
 * @file Vertex.h
 * @brief << brief description >>
 * @author clonker
 * @date 16.03.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <vector>
#include <readdy/common/macros.h>
#include <readdy/model/Particle.h>
#include <ostream>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(top)
NAMESPACE_BEGIN(graph)

/**
 * Struct representing a vertex in a topology-connectivity-graph
 */
struct Vertex {
    /**
     * edge in the graph (i.e., pointer to neighboring vertex)
     */
    using vertex_edge = Vertex *;

    /**
     * default constructor
     */
    Vertex() = default;

    /**
     * constructs a vertex to a graph
     * @param particleIndex the particle index this vertex belongs to
     * @param name named vertex, can be left empty and is then ignored
     * @param neighbors neighbors of the vertex (i.e., edges)
     */
    Vertex(std::size_t particleIndex, const std::string &name = "", const std::vector<vertex_edge> &neighbors = {})
            : neighbors(neighbors), particleIndex(particleIndex), name(name) {}

    /**
     * default destructor
     */
    virtual ~Vertex() = default;

    /**
     * the edges (i.e., pointers to neighboring vertices)
     */
    std::vector<vertex_edge> neighbors;
    /**
     * vertex' name, can be left empty and is then ignored
     */
    std::string name;
    /**
     * particle index in the topology this vertex belongs to
     */
    std::size_t particleIndex;

    bool operator==(const Vertex &rhs) const {
        return particleIndex == rhs.particleIndex;
    }

    friend std::ostream &operator<<(std::ostream &os, const Vertex &vertex) {
        os << "Vertex[name: " << vertex.name << " particleIndex: "
           << vertex.particleIndex << " neighbors=[";
        for (const auto neighbor : vertex.neighbors) {
            os << neighbor->particleIndex << ",";
        }
        os << "]]";
        return os;
    }

    bool operator!=(const Vertex &rhs) const {
        return !(rhs == *this);
    }
};

NAMESPACE_END(graph)
NAMESPACE_END(top)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
