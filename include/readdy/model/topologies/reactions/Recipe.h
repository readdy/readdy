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
    using graph_t = actions::TopologyReactionAction::graph_t;
    using vertex_ref = graph_t::vertex_ref;
    using vertex_cref = graph_t::vertex_cref;
    using label_edge = graph_t::label_edge;
    using label_vertex = graph_t::label;

    Recipe(GraphTopology& topology);

    Recipe(Recipe &&) = default;

    Recipe &operator=(Recipe &&) = default;

    Recipe(const Recipe &) = delete;

    Recipe &operator=(const Recipe &) = delete;

    Recipe &changeParticleType(const label_vertex &of, const particle_type_type &to);

private:
    GraphTopology &topology;
    reaction_operations steps;
};

NAMESPACE_END(reactions)
NAMESPACE_END(top)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
