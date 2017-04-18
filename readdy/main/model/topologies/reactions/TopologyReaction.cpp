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
 * @file TopologyReaction.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 03.04.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <readdy/model/topologies/reactions/TopologyReaction.h>

namespace readdy {
namespace model {
namespace top {
namespace reactions {

TopologyReaction::TopologyReaction(const reaction_function& reaction_function, const rate_function &rate_function)
        : reaction_function_(reaction_function), rate_function_(rate_function) { }

double TopologyReaction::rate(const GraphTopology &topology) const {
    return rate_function_(topology);
}

TopologyReaction::reaction_recipe TopologyReaction::operations(GraphTopology &topology) const {
    return reaction_function_(topology);
}

const bool TopologyReaction::raises_if_invalid() const {
    return mode_.flags.test(mode::raise_or_rollback_flag);
}

void TopologyReaction::raise_if_invalid() {
    mode_.raise();
}

const bool TopologyReaction::rolls_back_if_invalid() const {
    return !raises_if_invalid();
}

void TopologyReaction::roll_back_if_invalid() {
    mode_.rollback();
}

const bool TopologyReaction::expects_connected_after_reaction() const {
    return mode_.flags.test(mode::expect_connected_or_create_children_flag);
}

void TopologyReaction::expect_connected_after_reaction() {
    mode_.expect_connected();
}

const bool TopologyReaction::creates_child_topologies_after_reaction() const {
    return !expects_connected_after_reaction();
}

void TopologyReaction::create_child_topologies_after_reaction() {
    mode_.create_children();
}

TopologyReaction::TopologyReaction(const TopologyReaction::reaction_function &reaction_function, const double &rate)
        : TopologyReaction(reaction_function, [rate](const GraphTopology&) { return rate; }) {}

void Mode::raise() {
    flags[raise_or_rollback_flag] = true;
}

void Mode::rollback() {
    flags[raise_or_rollback_flag] = false;
}

void Mode::expect_connected() {
    flags[expect_connected_or_create_children_flag] = true;
}

void Mode::create_children() {
    flags[expect_connected_or_create_children_flag] = false;
}
}
}
}
}
