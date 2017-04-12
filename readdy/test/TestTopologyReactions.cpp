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
 * @file TestTopologyReactions.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 29.03.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <gtest/gtest.h>
#include <readdy/model/topologies/GraphTopology.h>

namespace {

TEST(TestTopologyReactions, ModeFlags) {
    using namespace readdy::model::top;
    reactions::TopologyReaction topologyReaction {[](const GraphTopology&) {
        reactions::TopologyReaction::reaction_recipe recipe;
        return recipe;
    }, [](const GraphTopology&) {
        return 0;
    }};
    topologyReaction.expect_connected_after_reaction();
    ASSERT_TRUE(topologyReaction.expects_connected_after_reaction());
    ASSERT_FALSE(topologyReaction.creates_child_topologies_after_reaction());

    topologyReaction.create_child_topologies_after_reaction();
    EXPECT_FALSE(topologyReaction.expects_connected_after_reaction());
    EXPECT_TRUE(topologyReaction.creates_child_topologies_after_reaction());

    topologyReaction.raise_if_invalid();
    EXPECT_TRUE(topologyReaction.raises_if_invalid());
    EXPECT_FALSE(topologyReaction.rolls_back_if_invalid());

    topologyReaction.roll_back_if_invalid();
    EXPECT_FALSE(topologyReaction.raises_if_invalid());
    EXPECT_TRUE(topologyReaction.rolls_back_if_invalid());
}

}
