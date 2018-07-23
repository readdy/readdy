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
 * @copyright GPL-3
 */
#pragma once

#include <readdy/common/macros.h>
#include "TopologyReactionAction.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(top)
NAMESPACE_BEGIN(reactions)
NAMESPACE_BEGIN(actions)

class TopologyReactionActionFactory {
public:
    /**
     * type of the graph topology's graph
     */
    using topology_graph = TopologyReactionAction::topology_graph;
    /**
     * reference to a topology reaction action
     */
    using action_ref = std::unique_ptr<TopologyReactionAction>;
    /**
     * reference of a vertex
     */
    using vertex = topology_graph::vertex_ref;
    /**
     * an edge in the graph
     */
    using edge = TopologyReactionAction::edge;

    /**
     * creates an action that changes the particle type of one of the particles contained in the topology
     * @param topology the topology
     * @param v the vertex to which the particle belongs whose type should be changed
     * @param type_to the target type
     * @return a unique pointer to the action
     */
    virtual action_ref createChangeParticleType(GraphTopology *topology, const vertex &v,
                                                const ParticleTypeId &type_to) const = 0;

    /**
     * creates an action that adds an edge in the topology
     * @param topology the topology
     * @param edge which edge to add
     * @return a unique pointer to the action
     */
    action_ref createAddEdge(GraphTopology *const topology, const edge &edge) const {
        return std::make_unique<AddEdge>(topology, edge);
    };

    /**
     * creates an action that removes an edge in the topology
     * @param topology the topology
     * @param edge which edge to remove
     * @return a unique pointer to the action
     */
    action_ref createRemoveEdge(GraphTopology *const topology, const edge &edge) const {
        return std::make_unique<RemoveEdge>(topology, edge);
    };

    /**
     * creates an action that changes the topology type
     * @param topology the topology
     * @param type_to the target type
     * @return a unique pointer to the action
     */
    virtual action_ref createChangeTopologyType(GraphTopology *topology, const std::string& type_to) const = 0;

};

NAMESPACE_END(actions)
NAMESPACE_END(reactions)
NAMESPACE_END(top)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
