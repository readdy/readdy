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
 * @file Graph.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 16.03.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <readdy/model/topologies/connectivity/Graph.h>
#include <algorithm>
#include <sstream>

namespace readdy {
namespace model {
namespace top {
namespace graph {

void Graph::addEdge(std::size_t v1, std::size_t v2) {
    if (v1 > vertices().size() || v2 > vertices().size() || v1 == v2) {
        throw std::invalid_argument(
                "When adding an edge the start and end vertices must not be equal and the "
                        "indices be contained in the vertices vector. v1=" + std::to_string(v1) + ", v2=" +
                std::to_string(v2) + ", #vertices=" + std::to_string(vertices().size()));
    }
    auto &vertex1 = vertices_.at(v1);
    auto &vertex2 = vertices_.at(v2);
    vertex1.neighbors.push_back(&vertex2);
    vertex2.neighbors.push_back(&vertex1);
}

void Graph::addEdge(std::size_t v1, const std::string &v2) {
    if (v1 > vertices().size()) {
        throw std::invalid_argument("the provided vertex index " + std::to_string(v1) +
                                    " was larger than the current number of vertices (" +
                                    std::to_string(vertices().size()) + ")!");
    }
    auto findIt = namedVertices.find(v2);
    if (findIt == namedVertices.end()) {
        throw std::invalid_argument("the provided vertex " + v2 + " is not contained in this graph.");
    }
    auto &vertex1 = vertices_.at(v1);
    auto vertex2Ptr = findIt->second;
    vertex1.neighbors.push_back(vertex2Ptr);
    vertex2Ptr->neighbors.push_back(&vertex1);
}

void Graph::addEdge(const std::string &v1, const std::string &v2) {
    auto findIt1 = namedVertices.find(v1);
    auto findIt2 = namedVertices.find(v2);
    if (findIt1 == namedVertices.end()) {
        throw std::invalid_argument("the provided vertex " + v1 + " is not contained in this graph.");
    }
    if (findIt2 == namedVertices.end()) {
        throw std::invalid_argument("the provided vertex " + v2 + " is not contained in this graph.");
    }
    findIt1->second->neighbors.push_back(findIt2->second);
    findIt2->second->neighbors.push_back(findIt1->second);
}

const std::vector<Vertex> &Graph::vertices() const {
    return vertices_;
}

void Graph::addVertex(Vertex &&v) {
    if (!v.name.empty()) {
        if (namedVertices.find(v.name) != namedVertices.end()) {
            throw std::invalid_argument("the named vertex \"" + v.name + "\" already existed!");
        }
        auto name = v.name;
        vertices_.push_back(std::move(v));
        namedVertices[name] = &vertices_.back();
    } else {
        vertices_.push_back(std::move(v));
    }
}

void Graph::addVertex(const Vertex &v) {
    addVertex(Vertex(v));
}

void Graph::removeVertex(std::size_t index) {
    if (index >= vertices_.size()) {
        throw std::invalid_argument("the index provided exceeded the number of vertices in this graph");
    }
    auto it = vertices_.begin() + index;
    removeNeighborsEdges(&*it);
    if(!it->name.empty()) {
        namedVertices.erase(namedVertices.find(it->name));
    }
    vertices_.erase(it);
}

void Graph::removeParticle(std::size_t particleIndex) {
    auto it = std::find_if(vertices_.begin(), vertices_.end(), [particleIndex](const Vertex &vertex) {
        return vertex.particleIndex == particleIndex;
    });
    if (it != vertices_.end()) {
        removeNeighborsEdges(&*it);
        if(!it->name.empty()) {
            namedVertices.erase(namedVertices.find(it->name));
        }
        vertices_.erase(it);
    } else {
        throw std::invalid_argument(
                "the vertex corresponding to the particle with topology index " + std::to_string(particleIndex) +
                " did not exist in the graph");
    }
}

void Graph::removeVertex(Vertex *vertex) {
    removeNeighborsEdges(vertex);
    auto it = std::find(vertices_.begin(), vertices_.end(), *vertex);
    if(!it->name.empty()) {
        namedVertices.erase(namedVertices.find(it->name));
    }
    if (it != vertices_.end()) {
        vertices_.erase(it);
    } else {
        std::stringstream ss;
        ss << *vertex;
        throw std::invalid_argument("the vertex " + ss.str() + " that was to be removed did not exist in the graph!");
    }
}

void Graph::removeVertex(const std::string &name) {
    decltype(namedVertices.begin()) it;
    if ((it = namedVertices.find(name)) == namedVertices.end()) {
        throw std::invalid_argument("the vertex \"" + name + "\" did not exist!");
    }
    removeVertex(it->second);
}

void Graph::removeNeighborsEdges(Vertex *vertex) {
    // need to remove myself from neighbors
    for (auto neighbor : vertex->neighbors) {
        auto it = std::find(neighbor->neighbors.begin(), neighbor->neighbors.end(), vertex);
        // I should always be my neighbors neighbor.
        assert(it != neighbor->neighbors.end());
        neighbor->neighbors.erase(it);
    }
}

void Graph::addEdge(const std::string &v1, std::size_t v2) {
    addEdge(v2, v1);
}

void Graph::removeEdge(Vertex *v1, Vertex *v2) {
    assert(v1 != v2);
    auto it1 = std::find(v1->neighbors.begin(), v1->neighbors.end(), v2);
    assert(it1 != v1->neighbors.end());
    auto it2 = std::find(v2->neighbors.begin(), v2->neighbors.end(), v1);
    assert(it2 != v2->neighbors.end());
    v1->neighbors.erase(it1);
    v2->neighbors.erase(it2);
}

void Graph::removeEdge(std::size_t v1, std::size_t v2) {
    assert(v1 < vertices_.size());
    assert(v2 < vertices_.size());
    removeEdge(&vertices_.at(v1), &vertices_.at(v2));
}

void Graph::removeEdge(std::size_t v1, const std::string &v2) {
    assert(v1 < vertices_.size());
    assert(namedVertices.find(v2) != namedVertices.end());
    removeEdge(&vertices_.at(v1), namedVertices.at(v2));
}

void Graph::removeEdge(const std::string &v1, std::size_t v2) {
    removeEdge(v2, v1);
}

void Graph::removeEdge(const std::string &v1, const std::string &v2) {
    assert(namedVertices.find(v1) != namedVertices.end());
    assert(namedVertices.find(v2) != namedVertices.end());
    removeEdge(namedVertices.at(v1), namedVertices.at(v2));
}


}
}
}
}