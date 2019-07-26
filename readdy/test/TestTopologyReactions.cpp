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
 * @copyright BSD-3
 */

#include <catch2/catch.hpp>

#include <readdy/model/topologies/GraphTopology.h>
#include <readdy/testing/KernelTest.h>
#include <readdy/testing/Utils.h>
#include <readdy/model/topologies/Utils.h>
#include <readdy/api/Simulation.h>

using namespace readdy;
using namespace readdytesting::kernel;

readdy::model::top::GraphTopology* setUpSmallTopology(readdy::model::Kernel* kernel, readdy::TopologyTypeId type) {
    auto &ctx = kernel->context();
    model::TopologyParticle x_0{0, 0, 0, ctx.particleTypes().idOf("Topology A")};
    model::TopologyParticle x_1{0, 0, 0, ctx.particleTypes().idOf("Topology A")};
    model::TopologyParticle x_2{0, 0, 0, ctx.particleTypes().idOf("Topology A")};
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

TEMPLATE_TEST_CASE("Test topology reactions.", "[topologies]", SingleCPU, CPU) {
    auto kernel = create<TestType>();
    auto &ctx = kernel->context();
    ctx.particleTypes().add("Topology A", 1.0, readdy::model::particleflavor::TOPOLOGY);
    ctx.particleTypes().add("Topology B", 1.0, readdy::model::particleflavor::TOPOLOGY);
    ctx.particleTypes().add("Topology Invalid Type", 1.0, readdy::model::particleflavor::TOPOLOGY);
    ctx.particleTypes().add("A", 1.0, readdy::model::particleflavor::NORMAL);

    ctx.topologyRegistry().configureBondPotential("Topology A", "Topology A", {10, 10});
    ctx.topologyRegistry().configureBondPotential("Topology A", "Topology B", {10, 10});
    ctx.topologyRegistry().configureBondPotential("Topology B", "Topology B", {10, 10});

    ctx.boxSize() = {{10, 10, 10}};

    SECTION("Mode flags") {
        using namespace readdy::model::top;

        reactions::StructuralTopologyReaction::reaction_function rfun = [](GraphTopology &t) -> reactions::Recipe {
            reactions::StructuralTopologyReaction::reaction_recipe recipe {t};
            return recipe;
        };
        reactions::StructuralTopologyReaction topologyReaction ("r", rfun, [](const GraphTopology &) {
            return 0;
        });
        topologyReaction.expect_connected_after_reaction();
        REQUIRE(topologyReaction.expects_connected_after_reaction());
        REQUIRE_FALSE(topologyReaction.creates_child_topologies_after_reaction());

        topologyReaction.create_child_topologies_after_reaction();
        REQUIRE_FALSE(topologyReaction.expects_connected_after_reaction());
        REQUIRE(topologyReaction.creates_child_topologies_after_reaction());

        topologyReaction.raise_if_invalid();
        REQUIRE(topologyReaction.raises_if_invalid());
        REQUIRE_FALSE(topologyReaction.rolls_back_if_invalid());

        topologyReaction.roll_back_if_invalid();
        REQUIRE_FALSE(topologyReaction.raises_if_invalid());
        REQUIRE(topologyReaction.rolls_back_if_invalid());
    }

    SECTION("Change a particle type") {
        auto tid = kernel->context().topologyRegistry().addType("Reactive Topology Type");
        auto topology = setUpSmallTopology(kernel.get(), tid);
        const auto &types = kernel->context().particleTypes();
        {
            auto reactionFunction = [&](model::top::GraphTopology &top) {
                model::top::reactions::Recipe recipe (top);
                recipe.changeParticleType(top.graph().vertices().begin(), types.idOf("Topology B"));
                return recipe;
            };
            model::top::reactions::StructuralTopologyReaction reaction {"r", reactionFunction, 5};
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
            model::top::reactions::StructuralTopologyReaction reaction {"r2", reactionFunction, rateFunction};
            reaction.expect_connected_after_reaction();
            reaction.raise_if_invalid();
            kernel->context().topologyRegistry().addStructuralReaction(tid, reaction);
        }
        const auto& reactions = kernel->context().topologyRegistry().structuralReactionsOf(tid);
        topology->updateReactionRates(reactions);
        REQUIRE(topology->rates().at(0) == 5); // Expected (constant) rate: 5
        REQUIRE(topology->rates().at(1) == 15); // Expected (function) rate: 15

        {
            auto result = reactions.at(0).execute(*topology, kernel.get());
            REQUIRE(result.empty()); // reaction is in-place, expect empty return vector
            auto particles = kernel->stateModel().getParticlesForTopology(*topology);
            auto v = topology->graph().vertices().begin();
            REQUIRE(particles[v->particleIndex].type() == types.idOf("Topology B"));
            // expect that the particle type in the graph representation and the particle data coincide
            REQUIRE(v->particleType() == particles[v->particleIndex].type());
        }
        {
            auto result = reactions.at(1).execute(*topology, kernel.get());
            REQUIRE(result.empty()); // reaction is in-place, expect empty return vector
            auto particles = kernel->stateModel().getParticlesForTopology(*topology);
            auto v = topology->graph().vertices().begin();
            REQUIRE(particles[v->particleIndex].type() == types.idOf("Topology A"));
            // expect that the particle type in the graph representation and the particle data coincide
            REQUIRE(v->particleType() == particles[v->particleIndex].type());
        }
    }
    SECTION("Add a named edge") {
        // add edge between first and last vertex, creating a circular structure x_2 -> x_0 <-> x_1 <-> x_2 <- x_0
        auto tid = kernel->context().topologyRegistry().addType("TA");
        auto topology = setUpSmallTopology(kernel.get(), tid);
        const auto &types = kernel->context().particleTypes();
        {
            auto reactionFunction = [&](model::top::GraphTopology &top) {
                model::top::reactions::Recipe recipe (top);
                recipe.addEdge(top.graph().vertices().begin(), std::prev(top.graph().vertices().end()));
                return recipe;
            };
            model::top::reactions::StructuralTopologyReaction reaction {"r", reactionFunction, 5};
            reaction.expect_connected_after_reaction();
            reaction.raise_if_invalid();
            kernel->context().topologyRegistry().addStructuralReaction(tid, reaction);
        }
        const auto &reactions = kernel->context().topologyRegistry().structuralReactionsOf(tid);
        topology->updateReactionRates(reactions);
        auto result = reactions.back().execute(*topology, kernel.get());
        REQUIRE(result.empty()); // reaction is in-place, expect empty return vector
        auto v1 = topology->graph().vertices().begin();
        auto v2 = std::prev(topology->graph().vertices().end());
        REQUIRE(topology->graph().containsEdge(std::make_tuple(v1, v2)));
    }
    SECTION("Add edge by iterator") {
        auto tid = kernel->context().topologyRegistry().addType("TA");
        auto topology = setUpSmallTopology(kernel.get(), tid);
        {
            auto reactionFunction = [&](model::top::GraphTopology &top) {
                model::top::reactions::Recipe recipe (top);
                recipe.addEdge(std::make_tuple(top.graph().vertices().begin(), --top.graph().vertices().end()));
                return recipe;
            };
            model::top::reactions::StructuralTopologyReaction reaction {"r", reactionFunction, 5};
            reaction.expect_connected_after_reaction();
            reaction.raise_if_invalid();
            kernel->context().topologyRegistry().addStructuralReaction(tid, reaction);
        }
        const auto &reactions = kernel->context().topologyRegistry().structuralReactionsOf(tid);
        topology->updateReactionRates(reactions);
        auto result = reactions.back().execute(*topology, kernel.get());
        REQUIRE(result.empty()); // reaction is in-place, expect empty return vector
        const auto &vertices = topology->graph().vertices();
        REQUIRE(topology->graph().containsEdge(std::make_tuple(vertices.begin(), std::prev(vertices.end()))));
    }
    SECTION("Remove edge straightforward") {
        kernel->context().topologyRegistry().addType("TA");
        auto topology = setUpSmallTopology(kernel.get(), kernel->context().topologyRegistry().idOf("TA"));

        {
            auto reactionFunction = [&](model::top::GraphTopology &top) {
                model::top::reactions::Recipe recipe (top);
                recipe.removeEdge(top.graph().vertices().begin(), ++top.graph().vertices().begin());
                return recipe;
            };
            model::top::reactions::StructuralTopologyReaction reaction {"r", reactionFunction, 5};
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
        REQUIRE(result.size() == 2);

        const auto& top1 = result.at(0);
        const auto& top2 = result.at(1);

        // first topology should have only 1 particle
        REQUIRE(top1.getNParticles() == 1);
        // second topology should have 2 particles
        REQUIRE(top2.getNParticles() == 2);

        SECTION("Second topology has one edge"){
            REQUIRE(top2.graph().vertices().begin()->neighbors().size() == 1);
            REQUIRE((++top2.graph().vertices().begin())->neighbors().size() == 1);
            REQUIRE(top2.graph().vertices().begin()->neighbors().at(0) == ++top2.graph().vertices().begin());
            REQUIRE((++top2.graph().vertices().begin())->neighbors().at(0) == top2.graph().vertices().begin());
        }
        SECTION("Particle indices are topology-relative and should begin with 0") {
            REQUIRE(top1.graph().vertices().begin()->particleIndex == 0);
            REQUIRE(top2.graph().vertices().begin()->particleIndex == 0);
            REQUIRE((++top2.graph().vertices().begin())->particleIndex == 1);
        }

        // check if particle mappings are still valid
        REQUIRE(top1.getParticles().at(0) == particles.at(0));
        REQUIRE(top2.getParticles().at(0) == particles.at(1));
        REQUIRE(top2.getParticles().at(1) == particles.at(2));
    }
    SECTION("Remove edge with rollback") {
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
            model::top::reactions::StructuralTopologyReaction reaction {"r", reactionFunction, 5};
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
        REQUIRE(result.empty());
        REQUIRE(graph.containsEdge(vertices.begin(), ++vertices.begin()));
        REQUIRE(graph.containsEdge(++vertices.begin(), --vertices.end()));
    }

    SECTION("Chain split up") {
        std::size_t n_chain_elements = 50;
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
            model::top::reactions::StructuralTopologyReaction reaction {"r", reactionFunction, rateFunction};
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
            REQUIRE(topPtr->getNParticles() == 1);
            REQUIRE(topPtr->graph().vertices().size() == 1);
        }
    }

    SECTION("Chain split up with decay") {
        std::size_t n_chain_elements = 50;
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
            model::top::reactions::StructuralTopologyReaction reaction {"r", reactionFunction, rateFunction};
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
            model::top::reactions::StructuralTopologyReaction reaction {"r2", reactionFunction, rateFunction};
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

        REQUIRE(kernel->stateModel().getTopologies().empty());
    }

    SECTION("Reaction types") {
        SECTION("TTFusion") {
            model::Context ctx;

            ctx.particleTypes().addTopologyType("X1", 0);
            ctx.particleTypes().addTopologyType("X2", 0);
            ctx.particleTypes().addTopologyType("Y", 0);
            ctx.particleTypes().addTopologyType("Z", 0);
            ctx.topologyRegistry().addType("T");
            ctx.topologyRegistry().addType("T2");

            ctx.topologyRegistry().configureBondPotential("Y", "Z", {0., .1});
            ctx.topologyRegistry().addSpatialReaction("connect: T(X1) + T(X2) -> T2(Y--Z)", 1e10, 1.);

            Simulation sim (std::move(kernel), ctx);

            std::string x1type = "X1"; 
            std::string x2type = "X2"; 
            std::string ttype = "T";
            REQUIRE(sim.context().topologyRegistry().spatialReactionsByType(x1type, ttype, x2type, ttype).size() == 1);
            REQUIRE(sim.context().topologyRegistry().spatialReactionsByType(x2type, ttype, x1type, ttype).size() == 1);

            auto p1 = sim.createTopologyParticle("X1", {0, 0, 0});
            auto p2 = sim.createTopologyParticle("X2", {0, 0, 0});
            sim.addTopology("T", {p1});
            sim.addTopology("T", {p2});

            auto topologies = sim.currentTopologies();

            REQUIRE(topologies.size() == 2);

            sim.run(1, 1e-3);

            topologies = sim.currentTopologies();

            REQUIRE(topologies.size() == 1);
            auto top = topologies.at(0);
            REQUIRE(top->type() == sim.context().topologyRegistry().idOf("T2"));
            REQUIRE(top->getNParticles() == 2);
            auto particles = top->fetchParticles();
            REQUIRE((particles[0].type() == sim.context().particleTypes().idOf("Y") ||
                        particles[0].type() == sim.context().particleTypes().idOf("Z")));
            if(particles[0].type() == sim.context().particleTypes().idOf("Y")) {
                REQUIRE(particles[1].type() == sim.context().particleTypes().idOf("Z"));
            } else {
                REQUIRE(particles[1].type() == sim.context().particleTypes().idOf("Y"));
            }
        }

        SECTION("TTSelfFusion") {
            model::Context ctx;

            ctx.particleTypes().addTopologyType("X", 0);
            ctx.particleTypes().addTopologyType("Y", 0);
            ctx.particleTypes().addTopologyType("Z", 0);
            ctx.topologyRegistry().addType("T");
            ctx.topologyRegistry().addType("T2");

            ctx.topologyRegistry().configureBondPotential("Y", "Z", {.0, .01});
            ctx.topologyRegistry().configureBondPotential("X", "X", {.0, .01});
            ctx.topologyRegistry().configureBondPotential("X", "Y", {.0, .01});
            ctx.topologyRegistry().configureBondPotential("X", "Z", {.0, .01});
            ctx.topologyRegistry().addSpatialReaction("connect: T(X) + T(X) -> T2(Y--Z) [self=true]", 1e10, 1.);

            Simulation sim (kernel->name(), ctx);

            auto p1 = sim.createTopologyParticle("X", {0, 0, -.01});
            auto p2 = sim.createTopologyParticle("X", {0, 0, 0});
            auto p3 = sim.createTopologyParticle("X", {0, 0, .01});
            auto t = sim.addTopology("T", {p1, p2, p3});
            t->graph().addEdgeBetweenParticles(0, 1);
            t->graph().addEdgeBetweenParticles(1, 2);

            sim.run(1, 1e-3);

            auto topologies = sim.currentTopologies();

            REQUIRE(topologies.size() == 1);
            auto top = topologies.at(0);
            REQUIRE(top->type() == sim.context().topologyRegistry().idOf("T2"));
            REQUIRE(top->getNParticles() == 3);
            auto particles = top->fetchParticles();
            for (auto ptype : {"X", "Y", "Z"}) {
                auto itX = std::find_if(particles.begin(), particles.end(), [&](const auto &p) {
                    return p.type() == sim.context().particleTypes().idOf(ptype);
                });
                INFO("at least one "+std::string(ptype)+" particle should be contained");
                REQUIRE(itX != particles.end());
                particles.erase(itX);
            }
        }
    }
}
