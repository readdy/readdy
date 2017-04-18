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
#include <readdy/testing/KernelTest.h>
#include <readdy/testing/Utils.h>

namespace {

using particle_t = readdy::model::Particle;
using topology_particle_t = readdy::model::TopologyParticle;

class TestTopologyReactions : public KernelTest {
public:
    readdy::model::top::GraphTopology *topology;
protected:
    virtual void SetUp() override {
        if (kernel->supportsTopologies()) {
            auto &ctx = kernel->getKernelContext();
            ctx.particle_types().add("Topology A", 1.0, 1.0, particle_t::FLAVOR_TOPOLOGY);
            ctx.particle_types().add("Topology B", 1.0, 1.0, particle_t::FLAVOR_TOPOLOGY);
            ctx.setBoxSize(10, 10, 10);
            topology_particle_t x_0{0, 0, 0, ctx.particle_types().id_of("Topology A")};
            topology_particle_t x_1{0, 0, 0, ctx.particle_types().id_of("Topology A")};
            topology_particle_t x_2{0, 0, 0, ctx.particle_types().id_of("Topology A")};
            topology = kernel->getKernelStateModel().addTopology({x_0, x_1, x_2});
            {
                auto it = topology->graph().vertices().begin();
                auto it2 = ++topology->graph().vertices().begin();
                topology->graph().addEdge(it, it2);
                std::advance(it, 1);
                std::advance(it2, 1);
                topology->graph().addEdge(it, it2);
            }
        }
    }

};

TEST(TestTopologyReactions, ModeFlags) {
    using namespace readdy;
    using namespace readdy::model::top;
    reactions::TopologyReaction topologyReaction{[](GraphTopology &t) {
        reactions::TopologyReaction::reaction_recipe recipe {t};
        return recipe;
    }, [](const GraphTopology &) {
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

TEST_P(TestTopologyReactions, ChangeParticleType) {
    using namespace readdy;
    if (!kernel->supportsTopologies()) {
        log::debug("kernel {} does not support topologies, thus skipping the test", kernel->getName());
        return;
    }
    topology->graph().setVertexLabel(topology->graph().vertices().begin(), "begin");
    {
        const auto &types = kernel->getKernelContext().particle_types();
        auto reactionFunction = [&](model::top::GraphTopology &top) {
            model::top::reactions::Recipe recipe (top);
            recipe.changeParticleType("begin", types.id_of("Topology B"));
            return recipe;
        };
        topology->addReaction({reactionFunction, 5});
    }
    {
        const auto &types = kernel->getKernelContext().particle_types();
        auto reactionFunction = [&](model::top::GraphTopology &top) {
            model::top::reactions::Recipe recipe (top);
            recipe.changeParticleType("begin", types.id_of("Topology B"));
            return recipe;
        };
        auto rateFunction = [&](const model::top::GraphTopology &top) {
            return 5;
        };
        topology->addReaction({reactionFunction, rateFunction});
    }
}

TEST_P(TestTopologyReactions, AddEdge) {
    using namespace readdy;
    if (!kernel->supportsTopologies()) {
        log::debug("kernel {} does not support topologies, thus skipping the test", kernel->getName());
        return;
    }
    // todo rollback/no rollback/fail/graceful
}

TEST_P(TestTopologyReactions, RemoveEdge) {
    using namespace readdy;
    if (!kernel->supportsTopologies()) {
        log::debug("kernel {} does not support topologies, thus skipping the test", kernel->getName());
        return;
    }
    // todo topology (decomposes/doesnt decompose) into child topologies (expected / not expected) with (rollback/norollback)
}

INSTANTIATE_TEST_CASE_P(TestTopologyReactions, TestTopologyReactions,
                        ::testing::ValuesIn(readdy::testing::getKernelsToTest()));

}
