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

#include <algorithm>
#include <sstream>
#include <unordered_set>
#include <readdy/model/topologies/graph/Graph.h>

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

const Graph::vertex_list &Graph::vertices() const {
    return vertices_;
}

void Graph::removeVertex(vertex_ref vertex) {
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
        if (!v->label.empty()) {
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

void Graph::removeNeighborsEdges(vertex_ref vertex) {
    for (auto neighbor : vertex->neighbors()) {
        neighbor->removeNeighbor(vertex);
    }
}

void Graph::removeEdge(vertex_ref v1, vertex_ref v2) {
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

void Graph::setVertexLabel(vertex_ref vertex, const std::string &label) {
    if(!label.empty()) {
        auto it = namedVertices.find(label);
        if (it == namedVertices.end()) {
            namedVertices[label] = vertex;
            vertex->label = label;
        } else {
            throw std::invalid_argument("the label " + label + " already existed in this topology!");
        }
    } else {
        auto it = namedVertices.begin();
        for(; it != namedVertices.end();) {
            if(it->second == vertex) {
                it = namedVertices.erase(it);
            } else {
                ++it;
            }
        }
    }
}

void Graph::addVertex(std::size_t particleIndex, particle_type_type particleType, const std::string &label) {
    if (!label.empty()) {
        if (namedVertices.find(label) != namedVertices.end()) {
            throw std::invalid_argument("the named vertex \"" + label + "\" already existed!");
        }
        vertices_.emplace_back(particleIndex, particleType, label);
        namedVertices[label] = --vertices().end();
    } else {
        vertices_.emplace_back(particleIndex, particleType);
    }
}

auto Graph::vertexItForParticleIndex(std::size_t particleIndex) -> decltype(vertices_.begin()) {
    auto it = std::find_if(vertices_.begin(), vertices_.end(), [particleIndex](const Vertex &vertex) {
        return vertex.particleIndex == particleIndex;
    });
    if (it != vertices_.end()) {
        return it;
    }
    return vertices_.end();
}

void Graph::addEdgeBetweenParticles(std::size_t particleIndex1, std::size_t particleIndex2) {
    auto it1 = vertexItForParticleIndex(particleIndex1);
    auto it2 = vertexItForParticleIndex(particleIndex2);
    if (it1 != vertices_.end() && it2 != vertices_.end()) {
        it1->addNeighbor(it2);
        it2->addNeighbor(it1);
    } else {
        throw std::invalid_argument("the particles indices did not exist...");
    }
}

void Graph::addEdge(vertex_ref v1, vertex_ref v2) {
    v1->addNeighbor(v2);
    v2->addNeighbor(v1);
}

Graph::vertex_list &Graph::vertices() {
    return vertices_;
}

Vertex &Graph::namedVertex(const std::string &name) {
    return *namedVertices.at(name);
}

Graph::vertex_ref Graph::firstVertex() {
    return vertices().begin();
}

Graph::vertex_ref Graph::lastVertex() {
    return --vertices().end();
}

bool Graph::isConnected() {
    std::for_each(vertices_.begin(), vertices_.end(), [](Vertex &v) { v.visited = false; });
    std::vector<vertex_ref> unvisited;
    unvisited.push_back(vertices_.begin());
    std::size_t n_visited = 0;
    while(!unvisited.empty()) {
        auto vertex = unvisited.back();
        unvisited.pop_back();
        if(!vertex->visited) {
            vertex->visited = true;
            ++n_visited;
            for (auto neighbor : vertex->neighbors()) {
                if (!neighbor->visited) {
                    unvisited.push_back(neighbor);
                }
            }
        }
    }
    return n_visited == vertices_.size();
}

void Graph::findNTuples(const edge_callback &tuple_callback,
                        const path_len_2_callback &triple_callback,
                        const path_len_3_callback &quadruple_callback) {
    for (auto &v : vertices_) {
        v.visited = false;
    }

    for (auto it = vertices_.begin(); it != vertices_.end(); ++it) {
        it->visited = true;
        auto v_type = it->particleType();
        auto v_idx = it->particleIndex;
        auto &neighbors = it->neighbors();
        for (auto it_neigh : neighbors) {
            auto vv_type = it_neigh->particleType();
            auto vv_idx = it_neigh->particleIndex;
            if (!it_neigh->visited) {
                log::debug("got type tuple ({}, {}) for particles {}, {}", v_type, vv_type, v_idx, vv_idx);
                // got edge (v, vv), now look for N(v)\{vv} and N(vv)\(N(v) + v)
                tuple_callback(std::tie(it, it_neigh));
                for (auto quad_it_1 : neighbors) {
                    // N(v)\{vv}
                    if (it_neigh != quad_it_1) {
                        auto vvv_type = quad_it_1->particleType();
                        auto vvv_idx = quad_it_1->particleIndex;
                        // got one end of the quadruple
                        for (auto quad_it_2 : it_neigh->neighbors()) {
                            // if this other neighbor is no neighbor of v and not v itself,
                            // we got the other end of the quadruple
                            auto no_circle =
                                    std::find(neighbors.begin(), neighbors.end(), quad_it_2) == neighbors.end();
                            if (quad_it_2 != it && no_circle) {
                                auto vvvv_type = quad_it_2->particleType();
                                auto vvvv_idx = quad_it_2->particleIndex;
                                log::debug("got type quadruple ({}, {}, {}, {}) for particles {}, {}, {}, {}", vvv_type, v_type,
                                           vv_type, vvvv_type, vvv_idx, v_idx, vv_idx, vvvv_idx);
                                quadruple_callback(std::tie(quad_it_1, it, it_neigh, quad_it_2));
                            }
                        }
                    }
                }
            }
            for (auto it_neigh2 : neighbors) {
                if (it_neigh2 != it_neigh && it_neigh->particleIndex < it_neigh2->particleIndex) {
                    auto vvv_type = it_neigh2->particleType();
                    auto vvv_idx = it_neigh2->particleIndex;
                    log::debug("got type triple ({}, {}, {}) for particles {}, {}, {}", vv_type, v_type, vvv_type,
                               vv_idx, v_idx, vvv_idx);
                    triple_callback(std::tie(it_neigh, it, it_neigh2));
                }
            }
        }
    }
}

std::tuple<std::vector<Graph::edge>, std::vector<Graph::path_len_2>, std::vector<Graph::path_len_3>>
Graph::findNTuples() {
    auto tuple = std::make_tuple(std::vector<edge>(), std::vector<path_len_2>(), std::vector<path_len_3>());
    findNTuples([&](const edge& tup) {
        std::get<0>(tuple).push_back(tup);
    }, [&](const path_len_2& path2) {
        std::get<1>(tuple).push_back(path2);
    }, [&](const path_len_3& path3) {
        std::get<2>(tuple).push_back(path3);
    });
    return tuple;
}

Graph::vertex_ref Graph::namedVertexPtr(const std::string &name) const {
    return namedVertices.at(name);
}

void Graph::addEdge(const edge &edge) {
    addEdge(std::get<0>(edge), std::get<1>(edge));
}

void Graph::removeEdge(const edge &edge) {
    removeEdge(std::get<0>(edge), std::get<1>(edge));
}

void Graph::addEdge(const Graph::label_edge &edge) {
    addEdge(std::get<0>(edge), std::get<1>(edge));
}

void Graph::removeEdge(const Graph::label_edge &edge) {
    removeEdge(std::get<0>(edge), std::get<1>(edge));
}

Graph::edge Graph::namedEdge(const Graph::label_edge &edge) const {
    return std::make_tuple(namedVertexPtr(std::get<0>(edge)), namedVertexPtr(std::get<1>(edge)));
}

}
}
}
}