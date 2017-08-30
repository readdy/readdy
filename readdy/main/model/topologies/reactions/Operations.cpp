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
 * @date 13.04.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <readdy/model/topologies/reactions/Operations.h>

namespace readdy {
namespace model {
namespace top {
namespace reactions {
namespace op {

Operation::action_ptr ChangeParticleType::create_action(topology_ref topology, factory_ref factory) const {
    return factory->createChangeParticleType(topology, _vertex, _type_to);
}

ChangeParticleType::ChangeParticleType(const Operation::vertex_ref &vertex, particle_type_type type_to)
        : _vertex(vertex), _type_to(type_to) {}

Operation::action_ptr AddEdge::create_action(topology_ref topology, Operation::factory_ref factory) const {
    return factory->createAddEdge(topology, _edge);
}

AddEdge::AddEdge(const Operation::edge &edge) : _edge(edge) {}

RemoveEdge::RemoveEdge(const Operation::edge &edge) : _edge(edge) {}

Operation::action_ptr RemoveEdge::create_action(topology_ref topology, factory_ref factory) const {
    return factory->createRemoveEdge(topology, _edge);
}

ChangeTopologyType::ChangeTopologyType(const std::string &type_to) : _type_to(type_to) {}

Operation::action_ptr
ChangeTopologyType::create_action(Operation::topology_ref topology, Operation::factory_ref factory) const {
    return factory->createChangeTopologyType(topology, _type_to);
}

}
}
}
}
}
