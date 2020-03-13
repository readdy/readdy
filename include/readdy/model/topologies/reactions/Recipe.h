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
 * @file TopologyReactionRecipeBuilder.h
 * @brief << brief description >>
 * @author clonker
 * @date 13.04.17
 * @copyright BSD-3
 */

#pragma once

#include <vector>

#include "Operations.h"
#include "TopologyReactionAction.h"

namespace readdy::model {
class Kernel;
namespace top::reactions {

class Recipe {
public:

    using reaction_operations = std::vector<op::Operation::OperationPtr>;
    using graph_topology = GraphTopology;

    explicit Recipe(graph_topology &topology) : _topology(topology) {};

    Recipe(Recipe &&) = default;

    Recipe &operator=(Recipe &&) = default;

    Recipe(const Recipe &) = default;

    Recipe &operator=(const Recipe &) = default;

    ~Recipe() = default;

    Recipe &changeParticleType(Graph::PersistentVertexIndex ref, const std::string &to);

    Recipe &changeParticleType(Graph::PersistentVertexIndex ref, const ParticleTypeId &to) {
        _steps.push_back(std::make_shared<op::ChangeParticleType>(ref, to));
        return *this;
    }

    Recipe &appendNewParticle(const std::vector<Graph::PersistentVertexIndex> &neighbors, const std::string &type,
                              const Vec3 &position);

    Recipe &changeParticlePosition(Graph::PersistentVertexIndex ref, const Vec3 &pos) {
        _steps.push_back(std::make_shared<op::ChangeParticlePosition>(ref, pos));
        return *this;
    }

    Recipe &addEdge(const Graph::Edge &edge) {
        _steps.push_back(std::make_shared<op::AddEdge>(edge));
        return *this;
    }

    Recipe &addEdge(Graph::PersistentVertexIndex v1, Graph::PersistentVertexIndex v2) {
        return addEdge(std::tie(v1, v2));
    }

    Recipe &removeEdge(const Graph::Edge &edge) {
        _steps.push_back(std::make_shared<op::RemoveEdge>(edge));
        return *this;
    }

    Recipe &removeEdge(Graph::PersistentVertexIndex v1, Graph::PersistentVertexIndex v2) {
        return removeEdge(std::tie(v1, v2));
    }

    Recipe &separateVertex(Graph::PersistentVertexIndex vertex);

    Recipe &changeTopologyType(const std::string &type) {
        _steps.push_back(std::make_shared<op::ChangeTopologyType>(type));
        return *this;
    }

    [[nodiscard]] const reaction_operations &steps() const {
        return _steps;
    }

    graph_topology &topology() {
        return _topology;
    }

    [[nodiscard]] const graph_topology &topology() const {
        return _topology;
    }

private:
    std::reference_wrapper<graph_topology> _topology;
    reaction_operations _steps;
};

}
}
