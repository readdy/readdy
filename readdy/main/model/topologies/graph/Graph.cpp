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
 * @copyright BSD-3
 */

#include <algorithm>
#include <sstream>
#include <unordered_set>
#include <readdy/model/topologies/graph/Graph.h>

namespace readdy::model::top::graph {

  bool Graph::areConnectedWithNOrLessEdges(std::size_t n, const Vertex &v1, const Vertex &v2) {
    if (n==0) return false;
    bool found = false;
    for (const auto neigh: v1.neighbors()) {
      if (neigh->particleIndex == v1.particleIndex) continue;
      if (neigh->particleIndex == v2.particleIndex) {
	return true;	
      }      
      found = areConnectedWithNOrLessEdges(n-1, *neigh, v2);
      if (found) return true;
    }
    return false;
  }
  
bool Graph::isConnected() const {
    std::for_each(_vertices.begin(), _vertices.end(), [](const Vertex &v) { v.visited = false; });

    auto &vert = const_cast<vertex_list&>(_vertices);

    std::vector<vertex_ref> unvisited;
    unvisited.emplace_back(vert.begin());
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

std::vector<Graph> Graph::connectedComponentsDestructive() {
    _dirty = true;
    std::vector<vertex_list> subVertexLists;
    {
        std::vector<std::vector<vertex_ref>> components;

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
