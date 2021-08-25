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
 * @file TopologyReactionRecipyBuilder.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 13.04.17
 * @copyright BSD-3
 */

#include <readdy/model/topologies/reactions/Recipe.h>
#include <readdy/model/Kernel.h>

namespace readdy::model::top::reactions {

Recipe &Recipe::changeParticleType(Graph::PersistentVertexIndex ref, const std::string &to) {
    return changeParticleType(ref, _topology.get().context().particleTypes().idOf(to));
}

Recipe &Recipe::separateVertex(Graph::PersistentVertexIndex vertex) {
    const auto &v = _topology.get().graph().vertices().at(vertex);
    std::for_each(v.neighbors().begin(), v.neighbors().end(), [this, vertex](auto neighbor) {
        this->removeEdge(std::make_tuple(vertex, neighbor));
    });
    return *this;
}

Recipe &Recipe::appendNewParticle(const std::vector<Graph::PersistentVertexIndex> &neighbors, const std::string &type,
                                  const Vec3 &position) {
    if(neighbors.empty()) throw std::invalid_argument("Cannot append new particle without specifying any neighbors.");
    auto typeId = _topology.get().context().particleTypes().idOf(type);
    _steps.push_back(std::make_shared<op::AppendParticle>(neighbors, typeId, position));
    return *this;
}

}