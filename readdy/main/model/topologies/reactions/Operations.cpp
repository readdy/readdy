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

#include <readdy/model/topologies/reactions/Operation.h>
#include <readdy/model/topologies/GraphTopology.h>

namespace readdy {
namespace model {
namespace top {
namespace reactions {
namespace op {

Operation::Operation(GraphTopology *const topology) : topology(topology){ }

Operation::graph_t::edge Operation::get_edge(const optional_edge &edge, const optional_label_edge &label_edge) const {
    // one needs to be true, the other one false
    assert(std::get<0>(edge) != std::get<0>(label_edge));
    if(std::get<0>(edge)) {
        return std::get<1>(edge);
    } else {
        return topology->graph().namedEdge(std::get<1>(label_edge));
    }
}

Operation::graph_t::vertex_ref Operation::get_vertex(const optional_vertex &v, const optional_label_vertex &lbl) const {
    // one needs to be true, the other one false
    assert(std::get<0>(v) != std::get<0>(lbl));
    if(std::get<0>(v)) {
        return std::get<1>(v);
    } else {
        return topology->graph().namedVertexPtr(std::get<1>(lbl));
    }
}

ChangeParticleType::ChangeParticleType(GraphTopology *const topology, const graph_t::vertex_ref &v,
                                       const particle_type_type &type_to)
        : Operation(topology), optional_vertex_({true, v}), type_to(type_to), previous_type(type_to){}

Operation::graph_t::vertex_ref ChangeParticleType::vertex() const {
    return get_vertex(optional_vertex_, optional_label_vertex_);
}

AddEdge::AddEdge(GraphTopology *const topology, const graph_t::edge& edge)
        : Operation(topology), optional_edge_({true, edge}) {}

void AddEdge::execute() {
    topology->graph().addEdge(get_edge(optional_edge_, optional_label_edge_));
}

void AddEdge::undo() {
    topology->graph().removeEdge(get_edge(optional_edge_, optional_label_edge_));
}

RemoveEdge::RemoveEdge(GraphTopology *const topology, const graph_t::edge& edge)
        : Operation(topology), optional_edge_({true, edge}) {}

void RemoveEdge::execute() {
    topology->graph().removeEdge(get_edge(optional_edge_, optional_label_edge_));
}

void RemoveEdge::undo() {
    topology->graph().addEdge(get_edge(optional_edge_, optional_label_edge_));
}



}
}
}
}
}
