/********************************************************************
 * Copyright © 2018 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * Redistribution and use in source and binary forms, with or       *
 * without modification, are permitted provided that the            *
 * following conditions are met:                                    *
 *  1. Redistributions of source code must retain the above         *
 *     copyright notice, this list of conditions and the            *
 *     following disclaimer.                                        *
 *  2. Redistributions in binary form must reproduce the above      *
 *     copyright notice, this list of conditions and the following  *
 *     disclaimer in the documentation and/or other materials       *
 *     provided with the distribution.                              *
 *  3. Neither the name of the copyright holder nor the names of    *
 *     its contributors may be used to endorse or promote products  *
 *     derived from this software without specific                  *
 *     prior written permission.                                    *
 *                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND           *
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,      *
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF         *
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE         *
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR            *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,     *
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,         *
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER *
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,      *
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)    *
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF      *
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
 ********************************************************************/


/**
 * << detailed description >>
 *
 * @file TestTopologyReactions.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 29.03.17
 * @copyright GPL-3
 */

#include <gtest/gtest.h>
#include <readdy/model/topologies/GraphTopology.h>
#include <readdy/testing/KernelTest.h>
#include <readdy/testing/Utils.h>
#include <readdy/model/topologies/Utils.h>
#include <readdy/api/Simulation.h>

namespace {

using particle_t = readdy::model::Particle;
using topology_particle_t = readdy::model::TopologyParticle;

class TestTopologyReactions : public KernelTest {
protected:
    void SetUp() override {
        if (kernel->supportsTopologies()) {
            auto &ctx = kernel->context();
            ctx.particleTypes().add("Topology A", 1.0, readdy::model::particleflavor::TOPOLOGY);
            ctx.particleTypes().add("Topology B", 1.0, readdy::model::particleflavor::TOPOLOGY);
            ctx.particleTypes().add("Topology Invalid Type", 1.0, readdy::model::particleflavor::TOPOLOGY);
            ctx.particleTypes().add("A", 1.0, readdy::model::particleflavor::NORMAL);

            ctx.topologyRegistry().configureBondPotential("Topology A", "Topology A", {10, 10});
            ctx.topologyRegistry().configureBondPotential("Topology A", "Topology B", {10, 10});
            ctx.topologyRegistry().configureBondPotential("Topology B", "Topology B", {10, 10});

            ctx.boxSize() = {{10, 10, 10}};
        }
    }

};

readdy::model::top::GraphTopology* setUpSmallTopology(readdy::model::Kernel* kernel, readdy::TopologyTypeId type) {
    auto &ctx = kernel->context();
    topology_particle_t x_0{0, 0, 0, ctx.particleTypes().idOf("Topology A")};
    topology_particle_t x_1{0, 0, 0, ctx.particleTypes().idOf("Topology A")};
    topology_particle_t x_2{0, 0, 0, ctx.particleTypes().idOf("Topology A")};
    auto topology = kernel->stateModel().addTopology(type, {x_0, x_1, x_2});
    {
        auto it = topology->graph().vertices().begin();
        auto it2 = ++topology->graph().vertices().begin();
        topology->graph().addEdge(it, it2);
        std::advance(it, 1);
        std::advance(it2, 1);
        topology->graph().addEdge(it, it2);
    }
    return topology;
}

TEST(TestTopologyReactions, ModeFlags) {
    using namespace readdy;
    using namespace readdy::model::top;

    reactions::StructuralTopologyReaction::reaction_function rfun = [](GraphTopology &t) -> reactions::Recipe {
        reactions::StructuralTopologyReaction::reaction_recipe recipe {t};
        return recipe;
    };
    reactions::StructuralTopologyReaction topologyReaction (rfun, [](const GraphTopology &) {
        return 0;
    });
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
        log::debug("kernel {} does not support topologies, thus skipping the test", kernel->name());
        return;
    }
    auto tid = kernel->context().topologyRegistry().addType("Reactive Topology Type");
    auto topology = setUpSmallTopology(kernel.get(), tid);
    const auto &types = kernel->context().particleTypes();
    {
        auto reactionFunction = [&](model::top::GraphTopology &top) {
            model::top::reactions::Recipe recipe (top);
            recipe.changeParticleType(top.graph().vertices().begin(), types.idOf("Topology B"));
            return recipe;
        };
        model::top::reactions::StructuralTopologyReaction reaction {reactionFunction, 5};
        reaction.expect_connected_after_reaction();
        reaction.raise_if_invalid();
        kernel->context().topologyRegistry().addStructuralReaction(tid, reaction);
    }
    {
        auto reactionFunction = [&](model::top::GraphTopology &top) {
            model::top::reactions::Recipe recipe (top);
            recipe.changeParticleType(top.graph().vertices().begin(), types.idOf("Topology A"));
            return recipe;
        };
        auto rateFunction = [&](const model::top::GraphTopology &top) {
            return 15;
        };
        model::top::reactions::StructuralTopologyReaction reaction {reactionFunction, rateFunction};
        reaction.expect_connected_after_reaction();
        reaction.raise_if_invalid();
        kernel->context().topologyRegistry().addStructuralReaction(tid, reaction);
    }
    const auto& reactions = kernel->context().topologyRegistry().structuralReactionsOf(tid);
    topology->updateReactionRates(reactions);
    EXPECT_EQ(topology->rates().at(0), 5) << "Expected (constant) rate: 5";
    EXPECT_EQ(topology->rates().at(1), 15) << "Expected (function) rate: 15";

    {
        auto result = reactions.at(0).execute(*topology, kernel.get());
        ASSERT_EQ(result.size(), 0) << "reaction is in-place, expect empty return vector";
        auto particles = kernel->stateModel().getParticlesForTopology(*topology);
        auto v = topology->graph().vertices().begin();
        ASSERT_EQ(particles[v->particleIndex].type(), types.idOf("Topology B"));
        ASSERT_EQ(v->particleType(), particles[v->particleIndex].type()) << "expect that the particle type in "
                            "the graph representation and the particle data coincide";
    }
    {
        auto result = reactions.at(1).execute(*topology, kernel.get());
        ASSERT_EQ(result.size(), 0) << "reaction is in-place, expect empty return vector";
        auto particles = kernel->stateModel().getParticlesForTopology(*topology);
        auto v = topology->graph().vertices().begin();
        ASSERT_EQ(particles[v->particleIndex].type(), types.idOf("Topology A"));
        ASSERT_EQ(v->particleType(), particles[v->particleIndex].type()) << "expect that the particle type in "
                            "the graph representation and the particle data coincide";
    }
}

TEST_P(TestTopologyReactions, GEXF) {
    using namespace readdy;
    if(!kernel->supportsTopologies()) {
        log::debug("kernel {} does not support topologies, thus skipping the test", kernel->name());
        return;
    }
    auto topology = setUpSmallTopology(kernel.get(), 0);
    auto middle = ++topology->graph().vertices().begin();
    auto gexf = model::top::util::to_gexf(topology->graph());
    EXPECT_TRUE(gexf.find("<node id=\"1\"") != std::string::npos) << "middle node should appear";
    EXPECT_TRUE(gexf.find("source=\"0\" target=\"1\"") != std::string::npos) << "first two vertices are connected";
}

TEST_P(TestTopologyReactions, AddEdgeNamed) {
    // add edge between first and last vertex, creating a circular structure x_2 -> x_0 <-> x_1 <-> x_2 <- x_0
    using namespace readdy;
    if (!kernel->supportsTopologies()) {
        log::debug("kernel {} does not support topologies, thus skipping the test", kernel->name());
        return;
    }
    auto tid = kernel->context().topologyRegistry().addType("TA");
    auto topology = setUpSmallTopology(kernel.get(), tid);
    const auto &types = kernel->context().particleTypes();
    {
        auto reactionFunction = [&](model::top::GraphTopology &top) {
            model::top::reactions::Recipe recipe (top);
            recipe.addEdge(top.graph().vertices().begin(), std::prev(top.graph().vertices().end()));
            return recipe;
        };
        model::top::reactions::StructuralTopologyReaction reaction {reactionFunction, 5};
        reaction.expect_connected_after_reaction();
        reaction.raise_if_invalid();
        kernel->context().topologyRegistry().addStructuralReaction(tid, reaction);
    }
    const auto &reactions = kernel->context().topologyRegistry().structuralReactionsOf(tid);
    topology->updateReactionRates(reactions);
    auto result = reactions.back().execute(*topology, kernel.get());
    ASSERT_EQ(result.size(), 0) << "reaction is in-place, expect empty return vector";
    EXPECT_TRUE(topology->graph().containsEdge(std::make_tuple(topology->graph().vertices().begin(),
                                                               std::prev(topology->graph().vertices().end()))));
}

TEST_P(TestTopologyReactions, AddEdgeIterator) {
    // add edge between first and last vertex, creating a circular structure x_2 -> x_0 <-> x_1 <-> x_2 <- x_0
    using namespace readdy;
    if (!kernel->supportsTopologies()) {
        log::debug("kernel {} does not support topologies, thus skipping the test", kernel->name());
        return;
    }
    auto tid = kernel->context().topologyRegistry().addType("TA");
    auto topology = setUpSmallTopology(kernel.get(), tid);
    {
        auto reactionFunction = [&](model::top::GraphTopology &top) {
            model::top::reactions::Recipe recipe (top);
            recipe.addEdge(std::make_tuple(top.graph().vertices().begin(), --top.graph().vertices().end()));
            return recipe;
        };
        model::top::reactions::StructuralTopologyReaction reaction {reactionFunction, 5};
        reaction.expect_connected_after_reaction();
        reaction.raise_if_invalid();
        kernel->context().topologyRegistry().addStructuralReaction(tid, reaction);
    }
    const auto &reactions = kernel->context().topologyRegistry().structuralReactionsOf(tid);
    topology->updateReactionRates(reactions);
    auto result = reactions.back().execute(*topology, kernel.get());
    ASSERT_EQ(result.size(), 0) << "reaction is in-place, expect empty return vector";
    const auto &vertices = topology->graph().vertices();
    EXPECT_TRUE(topology->graph().containsEdge(std::make_tuple(vertices.begin(), std::prev(vertices.end()))));
}


TEST_P(TestTopologyReactions, RemoveEdgeStraightforwardCase) {
    using namespace readdy;
    if (!kernel->supportsTopologies()) {
        log::debug("kernel {} does not support topologies, thus skipping the test", kernel->name());
        return;
    }
    kernel->context().topologyRegistry().addType("TA");
    auto topology = setUpSmallTopology(kernel.get(), kernel->context().topologyRegistry().idOf("TA"));

    {
        auto reactionFunction = [&](model::top::GraphTopology &top) {
            model::top::reactions::Recipe recipe (top);
            recipe.removeEdge(top.graph().vertices().begin(), ++top.graph().vertices().begin());
            return recipe;
        };
        model::top::reactions::StructuralTopologyReaction reaction {reactionFunction, 5};
        reaction.create_child_topologies_after_reaction();
        reaction.raise_if_invalid();
        kernel->context().topologyRegistry().addStructuralReaction("TA", reaction);
    }
    const auto &reactions = kernel->context().topologyRegistry().structuralReactionsOf("TA");
    topology->updateReactionRates(reactions);
    model::top::Topology::particle_indices particles;
    {
        std::copy(topology->getParticles().begin(), topology->getParticles().end(), std::back_inserter(particles));
    }
    auto result = reactions.back().execute(*topology, kernel.get());
    ASSERT_EQ(result.size(), 2);

    const auto& top1 = result.at(0);
    const auto& top2 = result.at(1);

    ASSERT_EQ(top1.getNParticles(), 1) << "first topology should have only 1 particle";
    ASSERT_EQ(top2.getNParticles(), 2) << "second topology should have 2 particles";
    ASSERT_EQ(top2.graph().vertices().begin()->neighbors().size(), 1) << "second topology has one edge";
    ASSERT_EQ((++top2.graph().vertices().begin())->neighbors().size(), 1) << "second topology has one edge";
    ASSERT_EQ(top2.graph().vertices().begin()->neighbors().at(0), ++top2.graph().vertices().begin())
                                << "second topology has one edge";
    ASSERT_EQ((++top2.graph().vertices().begin())->neighbors().at(0), top2.graph().vertices().begin())
                                << "second topology has one edge";
    ASSERT_EQ(top1.graph().vertices().begin()->particleIndex, 0)
                                << "particle indices are topology relative and should begin with 0";
    ASSERT_EQ(top2.graph().vertices().begin()->particleIndex, 0)
                                << "particle indices are topology relative and should begin with 0";
    ASSERT_EQ((++top2.graph().vertices().begin())->particleIndex, 1)
                                << "particle indices are topology relative and should begin with 0";

    // check if particle mappings are still valid
    ASSERT_EQ(top1.getParticles().at(0), particles.at(0));
    ASSERT_EQ(top2.getParticles().at(0), particles.at(1));
    ASSERT_EQ(top2.getParticles().at(1), particles.at(2));
}

TEST_P(TestTopologyReactions, RemoveEdgeRollback) {
    using namespace readdy;
    if (!kernel->supportsTopologies()) {
        log::debug("kernel {} does not support topologies, thus skipping the test", kernel->name());
        return;
    }

    auto &context = kernel->context();
    auto &toptypes = context.topologyRegistry();
    toptypes.addType("TA");
    auto topology = setUpSmallTopology(kernel.get(), toptypes.idOf("TA"));

    {
        auto reactionFunction = [&](model::top::GraphTopology &top) {
            model::top::reactions::Recipe recipe (top);
            recipe.removeEdge(top.graph().vertices().begin(), ++top.graph().vertices().begin());
            return recipe;
        };
        model::top::reactions::StructuralTopologyReaction reaction {reactionFunction, 5};
        reaction.expect_connected_after_reaction();
        reaction.roll_back_if_invalid();
        toptypes.addStructuralReaction("TA", reaction);
    }
    topology->updateReactionRates(toptypes.structuralReactionsOf("TA"));
    model::top::Topology::particle_indices particles;
    {
        std::copy(topology->getParticles().begin(), topology->getParticles().end(), std::back_inserter(particles));
    }
    std::vector<model::top::GraphTopology> result;
    {
        log::Level level (spdlog::level::err);
        result = toptypes.structuralReactionsOf("TA").back().execute(*topology, kernel.get());
    }
    const auto& graph = topology->graph();
    const auto& vertices = graph.vertices();
    EXPECT_TRUE(result.empty());
    EXPECT_TRUE(graph.containsEdge(vertices.begin(), ++vertices.begin()));
    EXPECT_TRUE(graph.containsEdge(++vertices.begin(), --vertices.end()));
}

TEST_P(TestTopologyReactions, SplitUpChain) {
    using namespace readdy;
    if (!kernel->supportsTopologies()) {
        log::debug("kernel {} does not support topologies, thus skipping the test", kernel->name());
        return;
    }

    std::size_t n_chain_elements = 50;
    auto &ctx = kernel->context();
    auto &toptypes = ctx.topologyRegistry();

    toptypes.addType("TA");

    ctx.boxSize() = {{10, 10, 10}};
    std::vector<readdy::model::TopologyParticle> topologyParticles;
    {
        topologyParticles.reserve(n_chain_elements);
        for (std::size_t i = 0; i < n_chain_elements; ++i) {
            const auto id = ctx.particleTypes().idOf("Topology A");
            topologyParticles.emplace_back(-5 + i * 10. / static_cast<readdy::scalar>(n_chain_elements), 0, 0, id);
        }
    }
    auto topology = kernel->stateModel().addTopology(toptypes.idOf("TA"), topologyParticles);
    {
        auto it = topology->graph().vertices().begin();
        auto it2 = ++topology->graph().vertices().begin();
        while(it2 != topology->graph().vertices().end()) {
            topology->graph().addEdge(it, it2);
            std::advance(it, 1);
            std::advance(it2, 1);
        }
    }

    {
        auto reactionFunction = [&](model::top::GraphTopology &top) {
            model::top::reactions::Recipe recipe (top);
            auto& vertices = top.graph().vertices();
            auto current_n_vertices = vertices.size();
            if(current_n_vertices > 1) {
                auto edge = readdy::model::rnd::uniform_int<>(0, static_cast<int>(current_n_vertices-2));
                auto it1 = vertices.begin();
                auto it2 = ++vertices.begin();
                for(int i = 0; i < edge; ++i) {
                    ++it1;
                    ++it2;
                }
                recipe.removeEdge(it1, it2);
            }
            return recipe;
        };
        auto rateFunction = [](const model::top::GraphTopology &top) {
            return top.getNParticles() > 1 ? top.getNParticles()/50. : 0;
        };
        model::top::reactions::StructuralTopologyReaction reaction {reactionFunction, rateFunction};
        reaction.create_child_topologies_after_reaction();
        reaction.roll_back_if_invalid();
        toptypes.addStructuralReaction("TA", reaction);
    }

    {
        auto integrator = kernel->actions().createIntegrator("EulerBDIntegrator", 1.0);
        auto forces = kernel->actions().calculateForces();
        auto topReactions = kernel->actions().evaluateTopologyReactions(1.0);

        std::size_t time = 0;
        std::size_t n_time_steps = 500;

        kernel->initialize();

        forces->perform();
        kernel->evaluateObservables(time);
        for(time = 1; time < n_time_steps; ++time) {
            integrator->perform();
            topReactions->perform();
            forces->perform();
            kernel->evaluateObservables(time);
        }
    }

    auto topologies = kernel->stateModel().getTopologies();
    for(const auto topPtr : topologies) {
        EXPECT_EQ(topPtr->getNParticles(), 1);
        EXPECT_EQ(topPtr->graph().vertices().size(), 1);
    }

}

TEST_P(TestTopologyReactions, SplitUpChainDecay) {
    using namespace readdy;
    if (!kernel->supportsTopologies()) {
        log::debug("kernel {} does not support topologies, thus skipping the test", kernel->name());
        return;
    }

    std::size_t n_chain_elements = 50;
    auto &ctx = kernel->context();
    auto &toptypes = ctx.topologyRegistry();
    toptypes.addType("TA");

    ctx.boxSize() = {{10, 10, 10}};
    std::vector<readdy::model::TopologyParticle> topologyParticles;
    {
        topologyParticles.reserve(n_chain_elements);
        for (std::size_t i = 0; i < n_chain_elements; ++i) {
            const auto id = ctx.particleTypes().idOf("Topology A");
            topologyParticles.emplace_back(-5 + i * 10. / static_cast<readdy::scalar>(n_chain_elements), 0, 0, id);
        }
    }
    auto topology = kernel->stateModel().addTopology(toptypes.idOf("TA"), topologyParticles);
    {
        auto it = topology->graph().vertices().begin();
        auto it2 = ++topology->graph().vertices().begin();
        while(it2 != topology->graph().vertices().end()) {
            topology->graph().addEdge(it, it2);
            std::advance(it, 1);
            std::advance(it2, 1);
        }
    }

    {
        // split reaction
        auto reactionFunction = [&](model::top::GraphTopology &top) {
            model::top::reactions::Recipe recipe (top);
            auto& vertices = top.graph().vertices();
            auto current_n_vertices = vertices.size();
            if(current_n_vertices > 1) {
                auto edge = readdy::model::rnd::uniform_int<>(0, static_cast<int>(current_n_vertices - 2));
                auto it1 = vertices.begin();
                auto it2 = ++vertices.begin();
                for(int i = 0; i < edge; ++i) {
                    ++it1;
                    ++it2;
                }
                recipe.removeEdge(it1, it2);
            }

            return recipe;
        };
        auto rateFunction = [](const model::top::GraphTopology &top) {
            return top.getNParticles() > 1 ? top.getNParticles()/50. : 0;
        };
        model::top::reactions::StructuralTopologyReaction reaction {reactionFunction, rateFunction};
        reaction.create_child_topologies_after_reaction();
        reaction.roll_back_if_invalid();

        toptypes.addStructuralReaction("TA", reaction);
    }
    {
        // decay reaction
        auto reactionFunction = [&](model::top::GraphTopology &top) {
            model::top::reactions::Recipe recipe (top);
            if(top.graph().vertices().size() == 1) {
                recipe.changeParticleType(top.graph().vertices().begin(),
                                          kernel->context().particleTypes().idOf("A"));
            } else {
                throw std::logic_error("this reaction should only be executed when there is exactly "
                                               "one particle in the topology");
            }
            return recipe;
        };
        auto rateFunction = [](const model::top::GraphTopology &top) {
            return top.getNParticles() > 1 ? 0 : 1;
        };
        model::top::reactions::StructuralTopologyReaction reaction {reactionFunction, rateFunction};
        reaction.create_child_topologies_after_reaction();
        reaction.roll_back_if_invalid();
        toptypes.addStructuralReaction("TA", reaction);
    }

    {
        auto integrator = kernel->actions().createIntegrator("EulerBDIntegrator", 1.0);
        auto forces = kernel->actions().calculateForces();
        auto topReactions = kernel->actions().evaluateTopologyReactions(1.0);

        std::size_t time = 0;
        std::size_t n_time_steps = 500;

        kernel->initialize();

        forces->perform();
        kernel->evaluateObservables(time);
        for(time = 1; time < n_time_steps; ++time) {
            integrator->perform();
            topReactions->perform();
            forces->perform();
            kernel->evaluateObservables(time);

        }
    }

    EXPECT_EQ(kernel->stateModel().getTopologies().size(), 0);
}

TEST_P(TestTopologyReactions, ChainDecayIntegrationTest) {
    using namespace readdy;
    if (!kernel->supportsTopologies()) {
        log::debug("kernel {} does not support topologies, thus skipping the test", kernel->name());
        return;
    }
    auto &ctx = kernel->context();
    ctx.boxSize() = {{150, 150, 150}};
    ctx.periodicBoundaryConditions() = {{true, true, true}};

    ctx.particleTypes().add("Decay", 1.);
    ctx.particleTypes().addTopologyType("T", 1.);
    ctx.particleTypes().add("whatever", 1.);
    ctx.reactions().add("decay: Decay ->", 1000.);
    ctx.reactions().add("whatever: whatever ->", 1000.);
    ctx.potentials().addBox("Decay", 10, {-70, -70, -70}, {130, 130, 130});
    ctx.potentials().addBox("T", 10, {-70, -70, -70}, {130, 130, 130});

    ctx.topologyRegistry().configureBondPotential("T", "T", {20, 2});
    ctx.topologyRegistry().addType("unstable");
    // decay reaction
    auto reactionFunction = [&](model::top::GraphTopology &top) {
        model::top::reactions::Recipe recipe (top);
        auto& vertices = top.graph().vertices();
        auto current_n_vertices = vertices.size();
        auto vIdx = readdy::model::rnd::uniform_int<>(0, static_cast<int>(current_n_vertices - 1));
        auto v = vertices.begin();
        std::advance(v, vIdx);
        recipe.separateVertex(v);
        recipe.changeParticleType(v, "Decay");
        return recipe;
    };
    auto rateFunction = [](const model::top::GraphTopology &top) {
        return .1;
    };
    model::top::reactions::StructuralTopologyReaction reaction {reactionFunction, rateFunction};
    ctx.topologyRegistry().addStructuralReaction("unstable", reaction);

    std::vector<readdy::model::TopologyParticle> topologyParticles;
    {
        topologyParticles.reserve(70);
        for (std::size_t i = 0; i < 70; ++i) {
            const auto id = ctx.particleTypes().idOf("T");
            topologyParticles.emplace_back(-5 + i * 10. / static_cast<readdy::scalar>(70), 0, 0, id);
        }
    }

    for(auto i = 0_z; i < 50; ++i) {
        kernel->stateModel().addParticle(model::Particle(0, 0, 0, ctx.particleTypes().idOf("whatever")));
    }

    auto topology = kernel->stateModel().addTopology(ctx.topologyRegistry().idOf("unstable"), topologyParticles);
    {
        auto it = topology->graph().vertices().begin();
        auto it2 = ++topology->graph().vertices().begin();
        while(it2 != topology->graph().vertices().end()) {
            topology->graph().addEdge(it, it2);
            std::advance(it, 1);
            std::advance(it2, 1);
        }
    }

    auto topsobs = kernel->observe().topologies(1);
    auto connection = kernel->connectObservable(topsobs.get());
    topsobs->callback() = [&](const model::observables::Topologies::result_type &result) {
        auto tops = kernel->stateModel().getTopologies();
        EXPECT_EQ(result.size(), tops.size());
        for(std::size_t i = 0; i < tops.size(); ++i) {
            auto top = tops.at(i);
            auto record = result.at(i);

            auto edges = top->graph().edges();
            EXPECT_EQ(edges.size(), record.edges.size());

            auto particles = top->fetchParticles();

            EXPECT_EQ(record.particleIndices.size(), particles.size());
            for(std::size_t j = 0; j < record.particleIndices.size(); ++j) {
                auto topParticles = top->getParticles();
                kernel->stateModel().toDenseParticleIndices(topParticles.begin(), topParticles.end());
                EXPECT_TRUE(std::find(topParticles.begin(), topParticles.end(),
                                      record.particleIndices.at(j)) != topParticles.end());
                EXPECT_TRUE(particles.at(j).type() == ctx.particleTypes().idOf("T"));
            }

            for(const auto &edge : record.edges) {
                auto it = std::find_if(edges.begin(), edges.end(), [&](const auto &e) {
                    auto p1Idx = std::get<0>(e)->particleIndex;
                    auto p2Idx = std::get<1>(e)->particleIndex;
                    return (p1Idx == std::get<0>(edge) && p2Idx == std::get<1>(edge))
                           || (p1Idx == std::get<1>(edge) && p2Idx == std::get<0>(edge));
                });
                EXPECT_TRUE(it != edges.end());
            }

        }
    };

    {
        auto integrator = kernel->actions().createIntegrator("EulerBDIntegrator", 1e-2);
        auto forces = kernel->actions().calculateForces();
        auto topReactions = kernel->actions().evaluateTopologyReactions(1e-2);
        auto reactions = kernel->actions().uncontrolledApproximation(1e-2);

        std::size_t time = 0;
        std::size_t n_time_steps = 10000;

        kernel->initialize();

        forces->perform();
        kernel->evaluateObservables(time);
        for(time = 1; time < n_time_steps; ++time) {
            integrator->perform();
            topReactions->perform();
            forces->perform();
            reactions->perform();
            kernel->evaluateObservables(time);
            for(auto topPtr : kernel->stateModel().getTopologies()) {
                // check that all topologies are just containing T particles and their edges are also fine
                for(const auto &p : topPtr->fetchParticles()) {
                    EXPECT_EQ(p.type(), ctx.particleTypes().idOf("T"));
                }
                for(auto edge : topPtr->graph().edges()) {
                    auto v1 = std::get<0>(edge);
                    auto v2 = std::get<1>(edge);
                    EXPECT_EQ(topPtr->particleForVertex(v1).type(), ctx.particleTypes().idOf("T"));
                    EXPECT_EQ(topPtr->particleForVertex(v2).type(), ctx.particleTypes().idOf("T"));
                }
            }
        }
    }
}

TEST_P(TestTopologyReactions, TTFusion) {
    using namespace readdy;
    if (!kernel->supportsTopologies()) {
        log::debug("kernel {} does not support topologies, thus skipping the test", kernel->name());
        return;
    }

    Simulation sim (kernel->name());

    sim.context().particleTypes().addTopologyType("X", 0);
    sim.context().particleTypes().addTopologyType("Y", 0);
    sim.context().particleTypes().addTopologyType("Z", 0);
    sim.context().topologyRegistry().addType("T");
    sim.context().topologyRegistry().addType("T2");

    sim.context().topologyRegistry().configureBondPotential("Y", "Z", {0., .1});
    sim.context().topologyRegistry().addSpatialReaction("connect: T(X) + T(X) -> T2(Y--Z) [self=true]", 1e10, 1.);

    auto p1 = sim.createTopologyParticle("X", {0, 0, 0});
    auto p2 = sim.createTopologyParticle("X", {0, 0, 0});
    sim.addTopology("T", {p1});
    sim.addTopology("T", {p2});

    sim.createLoop(1e-3).run(1);

    auto topologies = sim.currentTopologies();

    EXPECT_EQ(topologies.size(), 1);
    auto top = topologies.at(0);
    EXPECT_EQ(top->type(), sim.context().topologyRegistry().idOf("T2"));
    EXPECT_EQ(top->getNParticles(), 2);
    auto particles = top->fetchParticles();
    EXPECT_TRUE(particles[0].type() == sim.context().particleTypes().idOf("Y") ||
                particles[0].type() == sim.context().particleTypes().idOf("Z"));
    if(particles[0].type() == sim.context().particleTypes().idOf("Y")) {
        EXPECT_TRUE(particles[1].type() == sim.context().particleTypes().idOf("Z"));
    } else {
        EXPECT_TRUE(particles[1].type() == sim.context().particleTypes().idOf("Y"));
    }
}

TEST_P(TestTopologyReactions, TTSelfFusion) {
    using namespace readdy;
    if (!kernel->supportsTopologies()) {
        log::debug("kernel {} does not support topologies, thus skipping the test", kernel->name());
        return;
    }

    Simulation sim (kernel->name());

    sim.context().particleTypes().addTopologyType("X", 0);
    sim.context().particleTypes().addTopologyType("Y", 0);
    sim.context().particleTypes().addTopologyType("Z", 0);
    sim.context().topologyRegistry().addType("T");
    sim.context().topologyRegistry().addType("T2");

    sim.context().topologyRegistry().configureBondPotential("Y", "Z", {.0, .01});
    sim.context().topologyRegistry().configureBondPotential("X", "X", {.0, .01});
    sim.context().topologyRegistry().configureBondPotential("X", "Y", {.0, .01});
    sim.context().topologyRegistry().configureBondPotential("X", "Z", {.0, .01});
    sim.context().topologyRegistry().addSpatialReaction("connect: T(X) + T(X) -> T2(Y--Z) [self=true]", 1e10, 1.);

    auto p1 = sim.createTopologyParticle("X", {0, 0, -.01});
    auto p2 = sim.createTopologyParticle("X", {0, 0, 0});
    auto p3 = sim.createTopologyParticle("X", {0, 0, .01});
    auto t = sim.addTopology("T", {p1, p2, p3});
    t->graph().addEdgeBetweenParticles(0, 1);
    t->graph().addEdgeBetweenParticles(1, 2);

    sim.createLoop(1e-3).run(1);

    auto topologies = sim.currentTopologies();

    EXPECT_EQ(topologies.size(), 1);
    auto top = topologies.at(0);
    EXPECT_EQ(top->type(), sim.context().topologyRegistry().idOf("T2"));
    EXPECT_EQ(top->getNParticles(), 3);
    auto particles = top->fetchParticles();
    for (auto ptype : {"X", "Y", "Z"}) {
        auto itX = std::find_if(particles.begin(), particles.end(), [&](const auto &p) {
            return p.type() == sim.context().particleTypes().idOf(ptype);
        });
        ASSERT_NE(itX, particles.end()) << "at least one "+std::string(ptype)+" particle should be contained";
        particles.erase(itX);
    }
}

INSTANTIATE_TEST_CASE_P(TestTopologyReactionsKernelTests, TestTopologyReactions,
                        ::testing::ValuesIn(readdy::testing::getKernelsToTest()));
}
