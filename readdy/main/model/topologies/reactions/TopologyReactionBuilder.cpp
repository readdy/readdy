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
 * @file TopologyReactionBuilder.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 07.04.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <readdy/model/topologies/reactions/TopologyReactionBuilder.h>

namespace readdy {
namespace model {
namespace top {
namespace reactions {


TopologyReactionBuilder::TopologyReactionBuilder() {}

TopologyReaction TopologyReactionBuilder::build(const Mode &mode, const TopologyReaction::rate_function &rate) {
    auto operations_copy = operations;
    operations.clear();
    return TopologyReaction(
            [operations_copy](const GraphTopology &topology) { return operations_copy; }, rate, mode
    );
}


}
}
}
}
