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
    using graph_t = TopologyReactionAction::graph_t;
    using operation_ref = std::unique_ptr<TopologyReactionAction>;
    using vertex_t = graph_t::vertex_ref;
    using edge_t = TopologyReactionAction::edge;

    virtual operation_ref createChangeParticleType(GraphTopology *const topology, const vertex_t &v,
                                                  const particle_type_type &type_to) const = 0;

    operation_ref createAddEdge(GraphTopology *const topology, const edge_t &edge) const {
        return std::make_unique<AddEdge>(topology, edge);
    };

    operation_ref createRemoveEdge(GraphTopology *const topology, const edge_t &edge) const {
        return std::make_unique<RemoveEdge>(topology, edge);
    };

};

NAMESPACE_END(actions)
NAMESPACE_END(reactions)
NAMESPACE_END(top)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
