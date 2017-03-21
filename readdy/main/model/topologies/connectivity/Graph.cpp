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

void Graph::addEdge(const std::string &v1, const std::string &v2) {
    auto findIt1 = namedVertices.find(v1);
    auto findIt2 = namedVertices.find(v2);
    if (findIt1 == namedVertices.end()) {
        throw std::invalid_argument("the provided vertex " + v1 + " is not contained in this graph.");
    }
    if (findIt2 == namedVertices.end()) {
        throw std::invalid_argument("the provided vertex " + v2 + " is not contained in this graph.");
    }
    findIt1->second->addNeighbor(findIt2->second);
    findIt2->second->addNeighbor(findIt1->second);
}

const Graph::vertices_t &Graph::vertices() const {
    return vertices_;
}

void Graph::removeVertex(vertices_t::iterator vertex) {
    removeNeighborsEdges(vertex);
    if (!vertex->label.empty()) {
        namedVertices.erase(namedVertices.find(vertex->label));
    }
    vertices_.erase(vertex);
}

void Graph::removeParticle(std::size_t particleIndex) {
    auto v = vertexItForParticleIndex(particleIndex);
    if (v != vertices_.end()) {
        removeNeighborsEdges(v);
        if(!v->label.empty()) {
            namedVertices.erase(namedVertices.find(v->label));
        }
        vertices_.erase(v);
    } else {
        throw std::invalid_argument(
                "the vertex corresponding to the particle with topology index " + std::to_string(particleIndex) +
                " did not exist in the graph");
    }
}

void Graph::removeVertex(const std::string &name) {
    decltype(namedVertices.begin()) it;
    if ((it = namedVertices.find(name)) == namedVertices.end()) {
        throw std::invalid_argument("the vertex \"" + name + "\" did not exist!");
    }
    removeVertex(it->second);
}

void Graph::removeNeighborsEdges(vertices_t::iterator vertex) {
    for (auto neighbor : vertex->neighbors()) {
        neighbor->removeNeighbor(vertex);
    }
}

void Graph::removeEdge(vertices_t::iterator v1, vertices_t::iterator v2) {
    assert(v1 != v2);
    v1->removeNeighbor(v2);
    v2->removeNeighbor(v1);
}

void Graph::removeEdge(const std::string &v1, const std::string &v2) {
    assert(namedVertices.find(v1) != namedVertices.end());
    assert(namedVertices.find(v2) != namedVertices.end());
    removeEdge(namedVertices.at(v1), namedVertices.at(v2));
}

const Vertex &Graph::namedVertex(const std::string &name) const {
    decltype(namedVertices.begin()) it;
    if ((it = namedVertices.find(name)) == namedVertices.end()) {
        throw std::invalid_argument("the requested vertex " + name + " did not exist.");
    }
    return *it->second;
}

const Vertex &Graph::vertexForParticleIndex(std::size_t particleIndex) const {
    auto it = std::find_if(vertices_.begin(), vertices_.end(), [particleIndex](const Vertex &vertex) {
        return vertex.particleIndex == particleIndex;
    });
    if (it != vertices_.end()) {
        return *it;
    } else {
        throw std::invalid_argument("graph did not contain the particle index " + std::to_string(particleIndex));
    }
}

void Graph::setVertexLabel(vertices_t::iterator vertex, const std::string &label) {
    auto it = namedVertices.find(label);
    if (it == namedVertices.end()) {
        namedVertices[label] = vertex;
        vertex->label = label;
    } else {
        throw std::invalid_argument("the label " + label + " already existed in this topology!");
    }
}

void Graph::addVertex(std::size_t particleIndex, const std::string &label) {
    if (!label.empty()) {
        if (namedVertices.find(label) != namedVertices.end()) {
            throw std::invalid_argument("the named vertex \"" + label + "\" already existed!");
        }
        vertices_.emplace_back(particleIndex, vertices_.size(), label);
        namedVertices[label] = --vertices().end();
    } else {
        vertices_.emplace_back(particleIndex, vertices_.size());
    }
}

auto Graph::vertexItForParticleIndex(std::size_t particleIndex) -> decltype(vertices_.begin()) {
    auto it = std::find_if(vertices_.begin(), vertices_.end(), [particleIndex](const Vertex &vertex) {
        return vertex.particleIndex == particleIndex;
    });
    if(it != vertices_.end()) {
        return it;
    }
    return vertices_.end();
}

void Graph::addEdgeBetweenParticles(std::size_t particleIndex1, std::size_t particleIndex2) {
    auto it1 = vertexItForParticleIndex(particleIndex1);
    auto it2 = vertexItForParticleIndex(particleIndex2);
    if(it1 != vertices_.end() && it2 != vertices_.end()) {
        it1->addNeighbor(it2);
        it2->addNeighbor(it1);
    } else {
        throw std::invalid_argument("the particles indices did not exist...");
    }
}

void Graph::addEdge(vertices_t::iterator v1, vertices_t::iterator v2) {
    v1->addNeighbor(v2);
    v2->addNeighbor(v1);
}

Graph::vertices_t &Graph::vertices() {
    return vertices_;
}

Vertex &Graph::namedVertex(const std::string &name) {
    return *namedVertices.at(name);
}

std::list<readdy::model::top::graph::Vertex>::iterator Graph::firstVertex() {
    return vertices().begin();
}

std::list<readdy::model::top::graph::Vertex>::iterator Graph::lastVertex() {
    return --vertices().end();
}

}
}
}
}