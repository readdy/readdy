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
 * @file Operations.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 05.04.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <readdy/model/topologies/reactions/TopologyReactionAction.h>
#include <readdy/model/topologies/GraphTopology.h>

namespace readdy {
namespace model {
namespace top {
namespace reactions {
namespace actions {

TopologyReactionAction::TopologyReactionAction(GraphTopology *const topology) : topology(topology){ }

ChangeParticleType::ChangeParticleType(GraphTopology *const topology, const vertex &v,
                                       const ParticleTypeId &type_to)
        : TopologyReactionAction(topology), _vertex(v), type_to(type_to), previous_type(type_to){}

AddEdge::AddEdge(GraphTopology *const topology, const edge &edge)
        : TopologyReactionAction(topology), label_edge_(edge) {}

void AddEdge::execute() {
    topology->graph().addEdge(label_edge_);
}

void AddEdge::undo() {
    topology->graph().removeEdge(label_edge_);
}

RemoveEdge::RemoveEdge(GraphTopology *const topology, const edge& edge)
        : TopologyReactionAction(topology), label_edge_(edge) {}

void RemoveEdge::execute() {
    topology->graph().removeEdge(label_edge_);
}

void RemoveEdge::undo() {
    topology->graph().addEdge(label_edge_);
}


ChangeTopologyType::ChangeTopologyType(GraphTopology *topology, TopologyTypeId newType)
        : TopologyReactionAction(topology), _newType(newType){}

void ChangeTopologyType::execute() {
    _prevType = topology->type();
    topology->type() = _newType;
}

void ChangeTopologyType::undo() {
    topology->type() = _prevType;
}

}
}
}
}
}
