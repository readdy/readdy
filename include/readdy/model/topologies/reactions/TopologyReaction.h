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
 * @file TopologyReaction.h
 * @brief << brief description >>
 * @author clonker
 * @date 03.04.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <vector>
#include <memory>
#include <functional>

#include <readdy/common/macros.h>
#include <bitset>

#include "Operation.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(top)
class GraphTopology;
NAMESPACE_BEGIN(reactions)

/**
 * Struct holding the mode of the reaction:
 * - whether it should raise or roll back after a failed reaction
 * - whether it is expected to be still connected or if it should be "fissionated"
 */
struct Mode {
    static constexpr std::size_t raise_or_rollback_flag = 0;
    static constexpr std::size_t expect_connected_or_create_children_flag = 1;
    std::bitset<2> flags;
};

class TopologyReaction {
public:
    using mode = Mode;
    using reaction_operations = std::vector<op::Operation::OperationRef>;
    using reaction_recipe = reaction_operations;
    using reaction_function = std::function<reaction_recipe(const GraphTopology &)>;
    using rate_function = std::function<double(const GraphTopology &)>;

    TopologyReaction(const reaction_function &reaction_function, const rate_function &rate_function);

    double rate(const GraphTopology &topology) const;

    reaction_recipe operations(const GraphTopology &topology) const;

    const bool raises_if_invalid() const;

    void raise_if_invalid();

    const bool rolls_back_if_invalid() const;

    void roll_back_if_invalid();

    const bool expects_connected_after_reaction() const;

    void expect_connected_after_reaction();

    const bool creates_child_topologies_after_reaction() const;

    void create_child_topologies_after_reaction();

private:
    rate_function rate_function_;
    reaction_function reaction_function_;
    mode mode_;
};

NAMESPACE_END(reactions)
NAMESPACE_END(top)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
