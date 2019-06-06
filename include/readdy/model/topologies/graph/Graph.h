/********************************************************************
 * Copyright © 2018 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * Redistribution and use in source and binary forms, with or       *
 * without modification, are permitted provided that the            *
 * following conditions are met:                                    *
 *  1. Redistributions of source code must retain the above         *
 *     copyright notice, this list of conditions and the            *
 *     following disclaimer.                                        *
 *  2. Redistributions in binary form must reproduce the above      *
 *     copyright notice, this list of conditions and the following  *
 *     disclaimer in the documentation and/or other materials       *
 *     provided with the distribution.                              *
 *  3. Neither the name of the copyright holder nor the names of    *
 *     its contributors may be used to endorse or promote products  *
 *     derived from this software without specific                  *
 *     prior written permission.                                    *
 *                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND           *
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,      *
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF         *
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE         *
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR            *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,     *
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,         *
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER *
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,      *
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)    *
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF      *
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
 ********************************************************************/


/**
 * << detailed description >>
 *
 * @file Graph.h
 * @brief << brief description >>
 * @author clonker
 * @date 16.03.17
 * @copyright BSD-3
 */

#pragma once

#include <stdexcept>
#include <list>
#include <unordered_map>
#include "Vertex.h"

namespace readdy::model::top::graph {

class Graph {
public:

    using vertex = Vertex;
    using vertex_list = std::list<vertex>;
    using vertex_ref = vertex_list::iterator;
    using vertex_cref = vertex_list::const_iterator;

    using edge = std::tuple<vertex_ref, vertex_ref>;
    using cedge = std::tuple<vertex_cref, vertex_cref>;

    using path_len_2 = std::tuple<vertex_ref, vertex_ref, vertex_ref>;
    using cpath_len_2 = std::tuple<vertex_cref, vertex_cref, vertex_cref>;

    using path_len_3 = std::tuple<vertex_ref, vertex_ref, vertex_ref, vertex_ref>;
    using cpath_len_3 = std::tuple<vertex_cref, vertex_cref, vertex_cref, vertex_cref>;

    Graph() = default;

    explicit Graph(vertex_list vertexList) : _vertices(std::move(vertexList)) {}

    virtual ~Graph() = default;

    Graph(const Graph &) = delete;

    Graph &operator=(const Graph &) = delete;

    Graph(Graph &&) = default;

    Graph &operator=(Graph &&) = default;

    const vertex_list &vertices() const {
        return _vertices;
    }

    vertex_list &vertices() {
        return _vertices;
    }

    vertex_ref firstVertex() {
        return vertices().begin();
    }

    vertex_ref lastVertex() {
        return --vertices().end();
    }

    vertex_ref toRef(const vertex &v) {
        auto it = std::find(std::begin(_vertices), std::end(_vertices), v);
        if (it != std::end(_vertices)) {
            return {it};
        }
        throw std::invalid_argument(fmt::format(
                "Provided vertex {} was not part of the graph, no ref could be created!", v
        ));
    }

    vertex_cref toRef(const vertex &v) const {
        auto it = std::find(std::begin(_vertices), std::end(_vertices), v);
        if (it != std::end(_vertices)) {
            return {it};
        }
        throw std::invalid_argument(fmt::format(
                "Provided vertex {} was not part of the graph, no ref could be created!", v
        ));
    }

    bool containsEdge(const cedge &edge) const {
        const auto& [v1, v2] = edge;
        const auto &v1Neighbors = v1->neighbors();
        const auto &v2Neighbors = v2->neighbors();
        return std::find(v1Neighbors.begin(), v1Neighbors.end(), v2) != v1Neighbors.end()
               && std::find(v2Neighbors.begin(), v2Neighbors.end(), v1) != v2Neighbors.end();
    }

    bool containsEdge(vertex_cref v1, vertex_cref v2) const {
        return containsEdge(std::tie(v1, v2));
    }

    const vertex &vertexForParticleIndex(std::size_t particleIndex) const {
        auto it = std::find_if(_vertices.begin(), _vertices.end(), [particleIndex](const vertex &vertex) {
            return vertex.particleIndex == particleIndex;
        });
        if (it != _vertices.end()) {
            return *it;
        }
        throw std::invalid_argument("graph did not contain the particle index " + std::to_string(particleIndex));
    }

    void addVertex(std::size_t particleIndex, ParticleTypeId particleType) {
        _vertices.emplace_back(particleIndex, particleType);
    }

    static void addEdge(vertex_ref v1, vertex_ref v2) {
        v1->addNeighbor(v2);
        v2->addNeighbor(v1);
    }

    void addEdge(const edge &edge) {
        addEdge(std::get<0>(edge), std::get<1>(edge));
    }

    void addEdgeBetweenParticles(std::size_t particleIndex1, std::size_t particleIndex2) {
        auto it1 = vertexItForParticleIndex(particleIndex1);
        auto it2 = vertexItForParticleIndex(particleIndex2);
        if (it1 != _vertices.end() && it2 != _vertices.end()) {
            it1->addNeighbor(it2);
            it2->addNeighbor(it1);
        } else {
            throw std::invalid_argument("the particles indices did not exist...");
        }
    }

    void removeEdge(vertex_ref v1, vertex_ref v2) {
        assert(v1 != v2);
        v1->removeNeighbor(v2);
        v2->removeNeighbor(v1);
    }

    void removeEdge(const edge &edge) {
        removeEdge(std::get<0>(edge), std::get<1>(edge));
    }

    void removeVertex(vertex_ref vertex) {
        removeNeighborsEdges(vertex);
        _vertices.erase(vertex);
    }

    void removeParticle(std::size_t particleIndex) {
        auto v = vertexItForParticleIndex(particleIndex);
        if (v != _vertices.end()) {
            removeNeighborsEdges(v);
            _vertices.erase(v);
        } else {
            throw std::invalid_argument(
                    "the vertex corresponding to the particle with topology index " + std::to_string(particleIndex) +
                    " did not exist in the graph");
        }
    }

    bool isConnected();

    std::vector<edge> edges() const {
        if(_dirty) {
            _edges.clear();
            findEdges([this](const edge &tup) {
                _edges.push_back(tup);
            });
            _dirty = false;
        }
        return _edges;
    };

    bool hasEdge(const edge &edge) const {
        const auto &[v1, v2] = edge;
        const auto &e = edges();
        auto it1 = std::find(v1->neighbors().begin(), v1->neighbors().end(), v2);
        auto it2 = std::find(v2->neighbors().begin(), v2->neighbors().end(), v1);
        return it1 != v1->neighbors().end() && it2 != v2->neighbors().end();
    }

    template<typename TupleCallback>
    void findEdges(const TupleCallback &edgeCallback) const;

    template<typename TupleCallback, typename TripleCallback, typename QuadrupleCallback>
    void findNTuples(const TupleCallback &tuple_callback,
                     const TripleCallback &triple_callback,
                     const QuadrupleCallback &quadruple_callback) const;

    std::tuple<std::vector<edge>, std::vector<path_len_2>, std::vector<path_len_3>>
    findNTuples() {
        auto tuple = std::make_tuple(std::vector<edge>(), std::vector<path_len_2>(), std::vector<path_len_3>());
        findNTuples([&](const edge &tup) {
            std::get<0>(tuple).push_back(tup);
        }, [&](const path_len_2 &path2) {
            std::get<1>(tuple).push_back(path2);
        }, [&](const path_len_3 &path3) {
            std::get<2>(tuple).push_back(path3);
        });
        return tuple;
    };

    /**
     * Returns the connected components, invalidates this graph
     * @return connected components
     */
    std::vector<Graph> connectedComponentsDestructive();

private:
    vertex_list _vertices{};

    mutable bool _dirty {true};
    mutable std::vector<edge> _edges;

    void removeNeighborsEdges(vertex_ref vertex) {
        std::for_each(std::begin(vertex->neighbors()), std::end(vertex->neighbors()), [vertex](const auto neighbor) {
            neighbor->removeNeighbor(vertex);
        });
    }

    auto vertexItForParticleIndex(std::size_t particleIndex) -> decltype(_vertices.begin()) {
        return std::find_if(_vertices.begin(), _vertices.end(), [particleIndex](const vertex &vertex) {
            return vertex.particleIndex == particleIndex;
        });
    }
};

template<typename TupleCallback>
void Graph::findEdges(const TupleCallback &edgeCallback) const {
    for (auto &v : _vertices) {
        v.visited = false;
    }

    for(auto it = _vertices.begin(); it != _vertices.end(); ++it) {
        it->visited = true;
        auto &neighbors = it->neighbors();
        for (auto it_neigh : neighbors) {
            if(!it_neigh->visited) {
                edgeCallback(std::tie(it, it_neigh));
            }
        }
    }
}

template<typename TupleCallback, typename TripleCallback, typename QuadrupleCallback>
void Graph::findNTuples(const TupleCallback &tuple_callback,
                 const TripleCallback &triple_callback,
                 const QuadrupleCallback &quadruple_callback) const {
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
                log::trace("got type tuple ({}, {}) for particles {}, {}", v_type, vv_type, v_idx, vv_idx);
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
                                log::trace("got type quadruple ({}, {}, {}, {}) for particles {}, {}, {}, {}", vvv_type, v_type,
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
                    log::trace("got type triple ({}, {}, {}) for particles {}, {}, {}", vv_type, v_type, vvv_type,
                               vv_idx, v_idx, vvv_idx);
                    triple_callback(std::tie(it_neigh, it, it_neigh2));
                }
            }
        }
    }
}


}
