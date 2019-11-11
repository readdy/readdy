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
 * @file OperationFactory.h
 * @brief << brief description >>
 * @author clonker
 * @date 05.04.17
 * @copyright BSD-3
 */
#pragma once

#include "TopologyReactionAction.h"

namespace readdy::model::top::reactions::actions {

class TopologyReactionActionFactory {
public:
    /**
     * reference to a topology reaction action
     */
    using ActionPtr = std::unique_ptr<TopologyReactionAction>;

    /**
     * creates an action that changes the particle type of one of the particles contained in the topology
     * @param topology the topology
     * @param v the vertex to which the particle belongs whose type should be changed
     * @param type_to the target type
     * @return a unique pointer to the action
     */
    virtual ActionPtr createChangeParticleType(GraphTopology *topology, const Graph::VertexIndex &v,
                                               const ParticleTypeId &type_to) const = 0;

    virtual ActionPtr createChangeParticlePosition(GraphTopology *topology, const Graph::VertexIndex &v, Vec3 position) const = 0;

    virtual ActionPtr createAppendParticle(GraphTopology *topology, const std::vector<Graph::VertexIndex> &neighbors,
                                           ParticleTypeId type, const Vec3 &position) const = 0;

    /**
     * creates an action that adds an edge in the topology
     * @param topology the topology
     * @param edge which edge to add
     * @return a unique pointer to the action
     */
    ActionPtr createAddEdge(GraphTopology *const topology, const Graph::Edge &edge) const {
        return std::make_unique<AddEdge>(topology, edge);
    };

    /**
     * creates an action that removes an edge in the topology
     * @param topology the topology
     * @param edge which edge to remove
     * @return a unique pointer to the action
     */
    ActionPtr createRemoveEdge(GraphTopology *const topology, const Graph::Edge &edge) const {
        return std::make_unique<RemoveEdge>(topology, edge);
    };

    /**
     * creates an action that changes the topology type
     * @param topology the topology
     * @param type_to the target type
     * @return a unique pointer to the action
     */
    virtual ActionPtr createChangeTopologyType(GraphTopology *topology, const std::string &type_to) const = 0;

};

}
