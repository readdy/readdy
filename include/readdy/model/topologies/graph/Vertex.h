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

#include <ostream>
#include <list>
#include <utility>
#include <vector>
#include <algorithm>
#include <readdy/common/common.h>
#include <readdy/model/Particle.h>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(top)
NAMESPACE_BEGIN(graph)

class Graph;

/**
 * Struct representing a vertex in a topology-connectivity-graph
 */
class Vertex {
public:
    /**
     * edge in the graph (i.e., pointer to neighboring vertex)
     */
    using vertex_ptr = std::list<Vertex>::iterator;
    using vertex_cptr = std::list<Vertex>::const_iterator;

    /**
     * default constructor
     */
    Vertex() = default;

    /**
     * constructs a vertex to a graph
     * @param particleIndex the particle index this vertex belongs to
     */
    Vertex(std::size_t particleIndex, ParticleTypeId particleType)
            : particleIndex(particleIndex), particleType_(particleType) {}

    Vertex(const Vertex &) = delete;

    Vertex &operator=(const Vertex &) = delete;

    Vertex(Vertex &&) = default;

    Vertex &operator=(Vertex &&) = default;

    /**
     * default destructor
     */
    virtual ~Vertex() = default;

    /**
     * particle index in the topology this vertex belongs to
     */
    std::size_t particleIndex {0};

    bool operator==(const Vertex &rhs) const {
        return particleIndex == rhs.particleIndex;
    }

    friend std::ostream &operator<<(std::ostream &os, const Vertex &vertex) {
        os << "Vertex[particleIndex: " << vertex.particleIndex << ", neighbors=[";
        for (const auto neighbor : vertex.neighbors_) {
            os << neighbor->particleIndex << ",";
        }
        os << "]]";
        return os;
    }

    bool operator!=(const Vertex &rhs) const {
        return !(rhs == *this);
    }

    void addNeighbor(const vertex_ptr &edge) {
        if (std::find(neighbors_.begin(), neighbors_.end(), edge) == neighbors_.end()) {
            neighbors_.push_back(edge);
        } else {
            log::warn("tried to add an already existing edge ({} - {})", particleIndex, edge->particleIndex);
        }
    }

    void removeNeighbor(const vertex_ptr &edge) {
        decltype(neighbors_.begin()) it;
        if ((it = std::find(neighbors_.begin(), neighbors_.end(), edge)) != neighbors_.end()) {
            neighbors_.erase(it);
        } else {
            log::warn("tried to remove a non existing edge {} - {}", particleIndex, edge->particleIndex);
        }
    }

    const std::vector<vertex_ptr> &neighbors() const {
        return neighbors_;
    }

    const ParticleTypeId &particleType() const {
        return particleType_;
    }

    void setParticleType(ParticleTypeId type) {
        particleType_ = type;
    }

    /**
     * flag if this vertex has been visited (for BFS/DFS)
     */
    bool visited {false};

private:
    friend class readdy::model::top::graph::Graph;

    /**
     * the edges (i.e., pointers to neighboring vertices)
     */
    std::vector<vertex_ptr> neighbors_{};

    ParticleTypeId particleType_ {0};
};

NAMESPACE_END(graph)
NAMESPACE_END(top)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
