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
#include <readdy/model/topologies/Utils.h>

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
            ctx.particle_types().add("Topology Invalid Type", 1.0, 1.0, particle_t::FLAVOR_TOPOLOGY);
            ctx.particle_types().add("A", 1.0, 1.0, particle_t::FLAVOR_NORMAL);

            ctx.configureTopologyBondPotential("Topology A", "Topology A", {10, 10});
            ctx.configureTopologyBondPotential("Topology A", "Topology B", {10, 10});
            ctx.configureTopologyBondPotential("Topology B", "Topology B", {10, 10});

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
    const auto &types = kernel->getKernelContext().particle_types();
    topology->graph().setVertexLabel(topology->graph().vertices().begin(), "begin");
    {
        auto reactionFunction = [&](model::top::GraphTopology &top) {
            model::top::reactions::Recipe recipe (top);
            recipe.changeParticleType(top.graph().vertices().begin(), types.id_of("Topology B"));
            return recipe;
        };
        model::top::reactions::TopologyReaction reaction {reactionFunction, 5};
        reaction.expect_connected_after_reaction();
        reaction.raise_if_invalid();
        topology->addReaction(reaction);
    }
    {
        auto reactionFunction = [&](model::top::GraphTopology &top) {
            model::top::reactions::Recipe recipe (top);
            recipe.changeParticleType("begin", types.id_of("Topology A"));
            return recipe;
        };
        auto rateFunction = [&](const model::top::GraphTopology &top) {
            return 15;
        };
        model::top::reactions::TopologyReaction reaction {reactionFunction, rateFunction};
        reaction.expect_connected_after_reaction();
        reaction.raise_if_invalid();
        topology->addReaction(reaction);
    }
    topology->updateReactionRates();
    auto& r1 = *topology->registeredReactions().begin();
    auto& r2 = *(topology->registeredReactions().begin() + 1);
    EXPECT_EQ(std::get<1>(r1), 5) << "Expected (constant) rate: 5";
    EXPECT_EQ(std::get<1>(r2), 15) << "Expected (function) rate: 15";

    {
        std::get<0>(r1).execute(*topology, kernel.get());
        auto particles = kernel->getKernelStateModel().getParticlesForTopology(*topology);
        auto v = topology->graph().namedVertexPtr("begin");
        ASSERT_EQ(particles[v->particleIndex].getType(), types.id_of("Topology B"));
        ASSERT_EQ(v->particleType(), particles[v->particleIndex].getType()) << "expect that the particle type in "
                            "the graph representation and the particle data coincide";
    }
    {
        std::get<0>(r2).execute(*topology, kernel.get());
        auto particles = kernel->getKernelStateModel().getParticlesForTopology(*topology);
        auto v = topology->graph().namedVertexPtr("begin");
        ASSERT_EQ(particles[v->particleIndex].getType(), types.id_of("Topology A"));
        ASSERT_EQ(v->particleType(), particles[v->particleIndex].getType()) << "expect that the particle type in "
                            "the graph representation and the particle data coincide";
    }
}

TEST_P(TestTopologyReactions, GEXF) {
    using namespace readdy;
    if(!kernel->supportsTopologies()) {
        log::debug("kernel {} does not support topologies, thus skipping the test", kernel->getName());
        return;
    }
    auto middle = ++topology->graph().vertices().begin();
    topology->graph().setVertexLabel(middle, "middle");
    auto gexf = model::top::util::to_gexf(topology->graph());
    EXPECT_TRUE(gexf.find("<node id=\"1\" label=\"middle\"") != std::string::npos) << "middle node should be labelled";
    EXPECT_TRUE(gexf.find("source=\"0\" target=\"1\"") != std::string::npos) << "first two vertices are connected";
}

TEST_P(TestTopologyReactions, AddEdge) {
    // add edge between first and last vertex, creating a circular structure x_2 -> x_0 <-> x_1 <-> x_2 <- x_0
    using namespace readdy;
    if (!kernel->supportsTopologies()) {
        log::debug("kernel {} does not support topologies, thus skipping the test", kernel->getName());
        return;
    }
    const auto &types = kernel->getKernelContext().particle_types();
    topology->graph().setVertexLabel(topology->graph().vertices().begin(), "begin");
    topology->graph().setVertexLabel(--topology->graph().vertices().end(), "end");
    {
        auto reactionFunction = [&](model::top::GraphTopology &top) {
            model::top::reactions::Recipe recipe (top);
            recipe.addEdge("begin", "end");
            return recipe;
        };
        model::top::reactions::TopologyReaction reaction {reactionFunction, 5};
        reaction.expect_connected_after_reaction();
        reaction.raise_if_invalid();
        topology->addReaction(reaction);
    }
    topology->updateReactionRates();
    std::get<0>(topology->registeredReactions().back()).execute(*topology, kernel.get());
    EXPECT_TRUE(topology->graph().containsEdge(std::make_tuple("begin", "end")));
}

TEST_P(TestTopologyReactions, RemoveEdge) {
    using namespace readdy;
    if (!kernel->supportsTopologies()) {
        log::debug("kernel {} does not support topologies, thus skipping the test", kernel->getName());
        return;
    }
    // todo topology (decomposes/doesnt decompose) into child topologies (expected / not expected) with (rollback/norollback)
}

INSTANTIATE_TEST_CASE_P(TestTopologyReactionsKernelTests, TestTopologyReactions,
                        ::testing::ValuesIn(readdy::testing::getKernelsToTest()));

}
