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
 * @file OperationFactory.h
 * @brief << brief description >>
 * @author clonker
 * @date 05.04.17
 * @copyright GNU Lesser General Public License v3.0
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
    virtual action_ref createChangeParticleType(GraphTopology *const topology, const vertex &v,
                                                  const particle_type_type &type_to) const = 0;

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
    virtual action_ref createChangeTopologyType(GraphTopology *const topology, const std::string& type_to) const = 0;

};

NAMESPACE_END(actions)
NAMESPACE_END(reactions)
NAMESPACE_END(top)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
