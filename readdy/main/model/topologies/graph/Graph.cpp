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
 * @file Graph.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 16.03.17
 * @copyright GPL-3
 */

#include <algorithm>
#include <sstream>
#include <unordered_set>
#include <readdy/model/topologies/graph/Graph.h>

namespace readdy {
namespace model {
namespace top {
namespace graph {

bool Graph::isConnected() {
    std::for_each(_vertices.begin(), _vertices.end(), [](Vertex &v) { v.visited = false; });
    std::vector<vertex_ref> unvisited;
    unvisited.emplace_back(_vertices.begin());
    std::size_t n_visited = 0;
    while(!unvisited.empty()) {
        auto vertex = unvisited.back();
        unvisited.pop_back();
        if(!vertex->visited) {
            vertex->visited = true;
            ++n_visited;
            for (auto neighbor : vertex->neighbors()) {
                if (!neighbor->visited) {
                    unvisited.emplace_back(neighbor);
                }
            }
        }
    }
    return n_visited == _vertices.size();
}

void Graph::findEdges(const Graph::edge_callback &edgeCallback) {
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

std::vector<Graph> Graph::connectedComponentsDestructive() {
    std::vector<vertex_list> subVertexLists;
    {
        std::vector<std::vector<vertex_cref>> components;

        std::for_each(_vertices.begin(), _vertices.end(), [](Vertex &v) { v.visited = false; });

        for(auto it = _vertices.begin(); it != _vertices.end(); ++it) {
            if(!it->visited) {
                // got a new component
                components.emplace_back();
                subVertexLists.emplace_back();

                auto& component = components.back();

                std::vector<vertex_ref> unvisitedInComponent;
                unvisitedInComponent.emplace_back(it);
                while (!unvisitedInComponent.empty()) {
                    auto& vertex = unvisitedInComponent.back();
                    unvisitedInComponent.pop_back();
                    if (!vertex->visited) {
                        vertex->visited = true;
                        {
                            component.emplace_back(vertex);
                        }
                        for (auto neighbor : vertex->neighbors()) {
                            if (!neighbor->visited) {
                                unvisitedInComponent.emplace_back(neighbor);
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
                    it_subLists->splice(it_subLists->end(), _vertices, vertex_ref);
                }
            }
        }
    }

    std::vector<Graph> subGraphs;
    subGraphs.reserve(subVertexLists.size());
    {
        for (auto &subVertexList : subVertexLists) {
            subGraphs.emplace_back(std::move(subVertexList));
        }
    }
    return std::move(subGraphs);
}

}
}
}
}
