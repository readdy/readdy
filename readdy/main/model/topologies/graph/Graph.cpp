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
    auto findIt1 = _vertex_label_mapping.find(v1);
    auto findIt2 = _vertex_label_mapping.find(v2);
    if (findIt1 == _vertex_label_mapping.end()) {
        throw std::invalid_argument("the provided vertex " + v1 + " is not contained in this graph.");
    }
    if (findIt2 == _vertex_label_mapping.end()) {
        throw std::invalid_argument("the provided vertex " + v2 + " is not contained in this graph.");
    }
    findIt1->second->addNeighbor(findIt2->second.data());
    findIt2->second->addNeighbor(findIt1->second.data());
}

const Graph::vertex_list &Graph::vertices() const {
    return _vertices;
}

void Graph::removeVertex(vertex_ref vertex) {
    removeNeighborsEdges(vertex);
    if (!vertex->label().empty()) {
        _vertex_label_mapping.erase(_vertex_label_mapping.find(vertex->label()));
    }
    _vertices.erase(vertex.data());
}

void Graph::removeParticle(std::size_t particleIndex) {
    auto v = vertexItForParticleIndex(particleIndex);
    if (v != _vertices.end()) {
        removeNeighborsEdges(v);
        if (!v->label().empty()) {
            _vertex_label_mapping.erase(_vertex_label_mapping.find(v->label()));
        }
        _vertices.erase(v);
    } else {
        throw std::invalid_argument(
                "the vertex corresponding to the particle with topology index " + std::to_string(particleIndex) +
                " did not exist in the graph");
    }
}

void Graph::removeVertex(const std::string &name) {
    decltype(_vertex_label_mapping.begin()) it;
    if ((it = _vertex_label_mapping.find(name)) == _vertex_label_mapping.end()) {
        throw std::invalid_argument("the vertex \"" + name + "\" did not exist!");
    }
    removeVertex(it->second);
}

void Graph::removeNeighborsEdges(vertex_ref vertex) {
    for (auto neighbor : vertex->neighbors()) {
        neighbor->removeNeighbor(vertex.data());
    }
}

void Graph::removeEdge(vertex_ref v1, vertex_ref v2) {
    assert(v1 != v2);
    v1->removeNeighbor(v2.data());
    v2->removeNeighbor(v1.data());
}

void Graph::removeEdge(const std::string &v1, const std::string &v2) {
    assert(_vertex_label_mapping.find(v1) != _vertex_label_mapping.end());
    assert(_vertex_label_mapping.find(v2) != _vertex_label_mapping.end());
    removeEdge(_vertex_label_mapping.at(v1), _vertex_label_mapping.at(v2));
}

const Vertex &Graph::namedVertex(const std::string &name) const {
    decltype(_vertex_label_mapping.begin()) it;
    if ((it = _vertex_label_mapping.find(name)) == _vertex_label_mapping.end()) {
        throw std::invalid_argument("the requested vertex " + name + " did not exist.");
    }
    return *it->second;
}

const Vertex &Graph::vertexForParticleIndex(std::size_t particleIndex) const {
    auto it = std::find_if(_vertices.begin(), _vertices.end(), [particleIndex](const Vertex &vertex) {
        return vertex.particleIndex == particleIndex;
    });
    if (it != _vertices.end()) {
        return *it;
    } else {
        throw std::invalid_argument("graph did not contain the particle index " + std::to_string(particleIndex));
    }
}

void Graph::setVertexLabel(vertex_ref vertex, const std::string &label) {
    if(!label.empty()) {
        auto it = _vertex_label_mapping.find(label);
        if (it == _vertex_label_mapping.end()) {
            _vertex_label_mapping.emplace(std::make_pair(label, vertex));
            vertex->_label = label;
        } else {
            throw std::invalid_argument("the label " + label + " already existed in this topology!");
        }
    } else {
        auto it = _vertex_label_mapping.begin();
        for(; it != _vertex_label_mapping.end();) {
            if(it->second == vertex) {
                it = _vertex_label_mapping.erase(it);
            } else {
                ++it;
            }
        }
    }
}

void Graph::addVertex(std::size_t particleIndex, particle_type_type particleType, const std::string &label) {
    if (!label.empty()) {
        if (_vertex_label_mapping.find(label) != _vertex_label_mapping.end()) {
            throw std::invalid_argument("the named vertex \"" + label + "\" already existed!");
        }
        _vertices.emplace_back(particleIndex, particleType, label);
        _vertex_label_mapping[label] = --vertices().end();
    } else {
        _vertices.emplace_back(particleIndex, particleType);
    }
}

auto Graph::vertexItForParticleIndex(std::size_t particleIndex) -> decltype(_vertices.begin()) {
    auto it = std::find_if(_vertices.begin(), _vertices.end(), [particleIndex](const Vertex &vertex) {
        return vertex.particleIndex == particleIndex;
    });
    if (it != _vertices.end()) {
        return it;
    }
    return _vertices.end();
}

void Graph::addEdgeBetweenParticles(std::size_t particleIndex1, std::size_t particleIndex2) {
    auto it1 = vertexItForParticleIndex(particleIndex1);
    auto it2 = vertexItForParticleIndex(particleIndex2);
    if (it1 != _vertices.end() && it2 != _vertices.end()) {
        it1->addNeighbor(it2);
        it2->addNeighbor(it1);
    } else {
        throw std::invalid_argument("the particles indices did not exist...");
    }
}

void Graph::addEdge(vertex_ref v1, vertex_ref v2) {
    v1->addNeighbor(v2.data());
    v2->addNeighbor(v1.data());
}

Graph::vertex_list &Graph::vertices() {
    return _vertices;
}

Vertex &Graph::namedVertex(const std::string &name) {
    return *_vertex_label_mapping.at(name);
}

Graph::vertex_ref Graph::firstVertex() {
    return vertices().begin();
}

Graph::vertex_ref Graph::lastVertex() {
    return --vertices().end();
}

bool Graph::isConnected() {
    std::for_each(_vertices.begin(), _vertices.end(), [](Vertex &v) { v.visited = false; });
    std::vector<vertex_ref> unvisited;
    unvisited.push_back(_vertices.begin());
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
    return n_visited == _vertices.size();
}

void Graph::findNTuples(const edge_callback &tuple_callback,
                        const path_len_2_callback &triple_callback,
                        const path_len_3_callback &quadruple_callback) {
    for (auto &v : _vertices) {
        v.visited = false;
    }

    for (auto it = _vertices.begin(); it != _vertices.end(); ++it) {
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

Graph::vertex_cref Graph::namedVertexPtr(const std::string &name) const {
    return _vertex_label_mapping.at(name);
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

Graph::cedge Graph::namedEdge(const Graph::label_edge &edge) const {
    return std::make_tuple(namedVertexPtr(std::get<0>(edge)), namedVertexPtr(std::get<1>(edge)));
}

const Graph::vertex_label_mapping &Graph::vertexLabelMapping() const {
    return _vertex_label_mapping;
}

Graph::vertex_label_mapping &Graph::vertexLabelMapping() {
    return _vertex_label_mapping;
}

Graph::edge Graph::namedEdge(const Graph::label_edge &edge) {
    return std::make_tuple(namedVertexPtr(std::get<0>(edge)), namedVertexPtr(std::get<1>(edge)));
}

Graph::vertex_ref Graph::namedVertexPtr(const Graph::label &name) {
    return _vertex_label_mapping.at(name);
}

std::vector<Graph> Graph::connectedComponentsDestructive() {
    std::vector<Graph> subGraphs;
    std::vector<vertex_list> subVertexLists;
    std::vector<vertex_label_mapping> subVertexLabelMappings;

    {
        std::vector<std::vector<vertex_cref>> components;

        std::for_each(_vertices.begin(), _vertices.end(), [](Vertex &v) { v.visited = false; });

        for(auto it = _vertices.begin(); it != _vertices.end(); ++it) {
            if(!it->visited) {
                // got a new component
                components.emplace_back();
                subVertexLists.emplace_back();
                subVertexLabelMappings.emplace_back();

                auto& component = components.back();
                auto& mapping = subVertexLabelMappings.back();

                {
                    component.emplace_back(it);
                    if (!it->_label.empty()) mapping[it->_label] = it;
                }

                std::vector<vertex_ref> unvisitedInComponent;
                unvisitedInComponent.emplace_back(it);
                while (!unvisitedInComponent.empty()) {
                    auto& vertex = unvisitedInComponent.back();
                    unvisitedInComponent.pop_back();
                    if (!vertex->visited) {
                        vertex->visited = true;
                        {
                            component.emplace_back(vertex);
                            if (!vertex->_label.empty()) mapping[vertex->_label] = vertex;
                        }
                        for (auto neighbor : vertex->neighbors()) {
                            if (!neighbor->visited) {
                                unvisitedInComponent.push_back(neighbor);
                            }
                        }
                    }
                }
            }
        }

        {
            // transfer vertices
            auto it_components = components.begin();
            auto it_subLists = subVertexLists.begin();
            for(; it_components != components.end(); ++it_components, ++it_subLists) {
                for(const auto& vertex_ref : *it_components) {
                    it_subLists->splice(it_subLists->end(), _vertices, vertex_ref.data());
                }
            }
        }
    }
    subGraphs.reserve(subVertexLists.size());
    {
        auto it_mappings = subVertexLabelMappings.begin();
        auto it_subLists = subVertexLists.begin();
        for(; it_mappings != subVertexLabelMappings.end(); ++it_mappings, ++it_subLists) {
            subGraphs.emplace_back(std::move(*it_subLists), std::move(*it_mappings));
        }
    }
    return subGraphs;
}

Graph::Graph(vertex_list vertexList, vertex_label_mapping vertexLabelMapping)
        : _vertices(std::move(vertexList)), _vertex_label_mapping(std::move(vertexLabelMapping)){}

bool Graph::containsEdge(const Graph::label_edge &edge) const {
    const auto& v1 = namedVertexPtr(std::get<0>(edge));
    const auto& v2 = namedVertexPtr(std::get<1>(edge));
    return containsEdge(std::make_tuple(v1, v2));
}

bool Graph::containsEdge(const Graph::cedge &edge) const {
    const auto& v1 = std::get<0>(edge);
    const auto& v2 = std::get<1>(edge);
    const auto& v1Neighbors = v1->neighbors();
    const auto& v2Neighbors = v2->neighbors();
    return std::find(v1Neighbors.begin(), v1Neighbors.end(), v2.data()) != v1Neighbors.end()
           && std::find(v2Neighbors.begin(), v2Neighbors.end(), v1.data()) != v2Neighbors.end();
}

}
}
}
}