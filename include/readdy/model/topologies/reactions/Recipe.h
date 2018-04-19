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
 * @file TopologyReactionRecipeBuilder.h
 * @brief << brief description >>
 * @author clonker
 * @date 13.04.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <vector>

#include <readdy/common/macros.h>

#include "Operations.h"
#include "TopologyReactionAction.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
class Kernel;
NAMESPACE_BEGIN(top)
NAMESPACE_BEGIN(reactions)

class Recipe {
public:

    using reaction_operations = std::vector<op::Operation::Ref>;
    using topology_graph = actions::TopologyReactionAction::topology_graph;
    using vertex_ref = topology_graph::vertex_ref;
    using vertex_cref = topology_graph::vertex_cref;
    using edge = topology_graph::edge;
    using graph_topology = GraphTopology;

    explicit Recipe(graph_topology &topology) : _topology(topology) {};

    Recipe(Recipe &&) = default;

    Recipe &operator=(Recipe &&) = default;

    Recipe(const Recipe &) = default;

    Recipe &operator=(const Recipe &) = default;

    ~Recipe() = default;

    Recipe &changeParticleType(const vertex_ref &ref, const std::string &to);

    Recipe &changeParticleType(const vertex_ref &ref, const particle_type_type &to) {
        _steps.push_back(std::make_shared<op::ChangeParticleType>(ref, to));
        return *this;
    }

    Recipe &addEdge(const edge &edge) {
        _steps.push_back(std::make_shared<op::AddEdge>(edge));
        return *this;
    }

    Recipe &addEdge(vertex_ref v1, vertex_ref v2) {
        return addEdge(std::tie(v1, v2));
    }

    Recipe &removeEdge(const edge &edge) {
        _steps.push_back(std::make_shared<op::RemoveEdge>(edge));
        return *this;
    }

    Recipe &removeEdge(vertex_ref v1, vertex_ref v2) {
        return removeEdge(std::tie(v1, v2));
    }

    Recipe &separateVertex(const vertex_ref &vertex) {
        std::for_each(vertex->neighbors().begin(), vertex->neighbors().end(), [this, &vertex](const auto &neighbor) {
            this->removeEdge(std::make_tuple(vertex, neighbor));
        });
        return *this;
    }
    
    Recipe &changeTopologyType(const std::string &type) {
        _steps.push_back(std::make_shared<op::ChangeTopologyType>(type));
        return *this;
    }

    const reaction_operations &steps() const {
        return _steps;
    }

    graph_topology &topology() {
        return _topology;
    }

    const graph_topology &topology() const {
        return _topology;
    }

private:
    std::reference_wrapper<graph_topology> _topology;
    reaction_operations _steps;
};

NAMESPACE_END(reactions)
NAMESPACE_END(top)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
