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

#include <catch2/catch_template_test_macros.hpp>

#include <readdy/model/topologies/GraphTopology.h>
#include <readdy/testing/KernelTest.h>
#include <readdy/testing/Utils.h>
#include <readdy/api/Simulation.h>

using namespace readdy;
using namespace readdytesting::kernel;

readdy::model::top::GraphTopology *setUpSmallTopology(readdy::model::Kernel *kernel, readdy::TopologyTypeId type) {
    auto &ctx = kernel->context();
    model::Particle x_0{0, 0, 0, ctx.particleTypes().idOf("Topology A")};
    model::Particle x_1{0, 0, 0, ctx.particleTypes().idOf("Topology A")};
    model::Particle x_2{0, 0, 0, ctx.particleTypes().idOf("Topology A")};
    auto topology = kernel->stateModel().addTopology(type, {x_0, x_1, x_2});
    {
        auto it = topology->graph().vertices().begin();
        auto it2 = ++topology->graph().vertices().begin();
        topology->addEdge({0}, {1});
        topology->addEdge({1}, {2});
    }
    return topology;
}

bool isNotConnected(readdy::model::top::Graph::iterator v) {
    return v->neighbors().empty();
}

readdy::scalar constant_rate_function(
        const readdy::model::top::GraphTopology &top1,
        const readdy::model::top::GraphTopology &top2
) {
    return 1e10;
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
            reactions::StructuralTopologyReaction::reaction_recipe recipe{t};
            return recipe;
        };
        reactions::StructuralTopologyReaction topologyReaction("r", rfun, [](const GraphTopology &) {
            return 0;
        });
        topologyReaction.expect_connected_after_reaction();
        REQUIRE(topologyReaction.expects_connected_after_reaction());
        REQUIRE_FALSE(topologyReaction.creates_child_topologies_after_reaction());

        topologyReaction.create_child_topologies_after_reaction();
        REQUIRE_FALSE(topologyReaction.expects_connected_after_reaction());
        REQUIRE(topologyReaction.creates_child_topologies_after_reaction());

    }

    SECTION("Change a particle type") {
        auto tid = kernel->context().topologyRegistry().addType("Reactive Topology Type");
        auto topology = setUpSmallTopology(kernel.get(), tid);
        const auto &types = kernel->context().particleTypes();
        {
            auto reactionFunction = [&](model::top::GraphTopology &top) {
                model::top::reactions::Recipe recipe(top);
                recipe.changeParticleType(top.graph().vertices().begin().persistent_index(), types.idOf("Topology B"));
                return recipe;
            };
            model::top::reactions::StructuralTopologyReaction reaction{"r", reactionFunction, 5};
            reaction.expect_connected_after_reaction();
            kernel->context().topologyRegistry().addStructuralReaction(tid, reaction);
        }
        {
            auto reactionFunction = [&](model::top::GraphTopology &top) {
                model::top::reactions::Recipe recipe(top);
                recipe.changeParticleType(top.graph().vertices().begin().persistent_index(), types.idOf("Topology A"));
                return recipe;
            };
            auto rateFunction = [&](const model::top::GraphTopology &top) {
                return 15;
            };
            model::top::reactions::StructuralTopologyReaction reaction{"r2", reactionFunction, rateFunction};
            reaction.expect_connected_after_reaction();
            kernel->context().topologyRegistry().addStructuralReaction(tid, reaction);
        }
        const auto &reactions = kernel->context().topologyRegistry().structuralReactionsOf(tid);
        topology->updateReactionRates(reactions);
        REQUIRE(topology->rates().at(0) == 5); // Expected (constant) rate: 5
        REQUIRE(topology->rates().at(1) == 15); // Expected (function) rate: 15

        {
            auto result = reactions.at(0).execute(*topology, kernel.get());
            REQUIRE(result.empty()); // reaction is in-place, expect empty return vector
            auto particles = kernel->stateModel().getParticlesForTopology(*topology);
            auto v = topology->graph().vertices().begin();
            REQUIRE(particles[v->data().particleIndex].type() == types.idOf("Topology B"));
        }
        {
            auto result = reactions.at(1).execute(*topology, kernel.get());
            REQUIRE(result.empty()); // reaction is in-place, expect empty return vector
            auto particles = kernel->stateModel().getParticlesForTopology(*topology);
            auto v = topology->graph().vertices().begin();
            REQUIRE(particles[v->data().particleIndex].type() == types.idOf("Topology A"));
        }
    }
    SECTION("Add a named edge") {
        // add edge between first and last vertex, creating a circular structure x_2 -> x_0 <-> x_1 <-> x_2 <- x_0
        auto tid = kernel->context().topologyRegistry().addType("TA");
        auto topology = setUpSmallTopology(kernel.get(), tid);
        const auto &types = kernel->context().particleTypes();
        {
            auto reactionFunction = [&](model::top::GraphTopology &top) {
                model::top::reactions::Recipe recipe(top);
                recipe.addEdge(top.graph().vertices().begin().persistent_index(), std::prev(top.graph().vertices().end()).persistent_index());
                return recipe;
            };
            model::top::reactions::StructuralTopologyReaction reaction{"r", reactionFunction, 5};
            reaction.expect_connected_after_reaction();
            kernel->context().topologyRegistry().addStructuralReaction(tid, reaction);
        }
        const auto &reactions = kernel->context().topologyRegistry().structuralReactionsOf(tid);
        topology->updateReactionRates(reactions);
        auto result = reactions.back().execute(*topology, kernel.get());
        REQUIRE(result.empty()); // reaction is in-place, expect empty return vector
        auto v1 = topology->graph().vertices().begin();
        auto v2 = std::prev(topology->graph().vertices().end());
        REQUIRE(topology->containsEdge(v1.persistent_index(), v2.persistent_index()));
    }
    SECTION("Add edge by iterator") {
        auto tid = kernel->context().topologyRegistry().addType("TA");
        auto topology = setUpSmallTopology(kernel.get(), tid);
        {
            auto reactionFunction = [&](model::top::GraphTopology &top) {
                model::top::reactions::Recipe recipe(top);
                recipe.addEdge(top.graph().vertices().begin().persistent_index(), (--top.graph().vertices().end()).persistent_index());
                return recipe;
            };
            model::top::reactions::StructuralTopologyReaction reaction{"r", reactionFunction, 5};
            reaction.expect_connected_after_reaction();
            kernel->context().topologyRegistry().addStructuralReaction(tid, reaction);
        }
        const auto &reactions = kernel->context().topologyRegistry().structuralReactionsOf(tid);
        topology->updateReactionRates(reactions);
        auto result = reactions.back().execute(*topology, kernel.get());
        REQUIRE(result.empty()); // reaction is in-place, expect empty return vector
        auto &vertices = topology->graph().vertices();
        REQUIRE(topology->graph().containsEdge(vertices.begin().persistent_index(), std::prev(vertices.end()).persistent_index()));
    }
    SECTION("Remove edge straightforward") {
        kernel->context().topologyRegistry().addType("TA");
        auto topology = setUpSmallTopology(kernel.get(), kernel->context().topologyRegistry().idOf("TA"));

        {
            auto reactionFunction = [&](model::top::GraphTopology &top) {
                model::top::reactions::Recipe recipe(top);
                recipe.removeEdge(top.graph().vertices().begin().persistent_index(), (++top.graph().vertices().begin()).persistent_index());
                return recipe;
            };
            model::top::reactions::StructuralTopologyReaction reaction{"r", reactionFunction, 5};
            reaction.create_child_topologies_after_reaction();
            kernel->context().topologyRegistry().addStructuralReaction("TA", reaction);
        }
        const auto &reactions = kernel->context().topologyRegistry().structuralReactionsOf("TA");
        topology->updateReactionRates(reactions);
        std::vector<model::Particle> particles = topology->fetchParticles();
        auto result = reactions.back().execute(*topology, kernel.get());
        REQUIRE(result.size() == 2);

        const auto &top1 = result.at(0);
        const auto &top2 = result.at(1);

        // first topology should have only 1 particle
        REQUIRE(top1.nParticles() == 1);
        // second topology should have 2 particles
        REQUIRE(top2.nParticles() == 2);

        SECTION("Second topology has one edge") {
            REQUIRE(top2.graph().vertices().begin()->neighbors().size() == 1);
            REQUIRE((++top2.graph().vertices().begin())->neighbors().size() == 1);
            REQUIRE(top2.graph().vertices().begin()->neighbors().at(0) == (++top2.graph().vertices().begin()).persistent_index());
            REQUIRE((++top2.graph().vertices().begin())->neighbors().at(0) == top2.graph().vertices().begin().persistent_index());
        }

        // check if particle mappings are still valid
        REQUIRE(top1.fetchParticles().at(0) == particles.at(0));
        REQUIRE(top2.fetchParticles().at(0) == particles.at(1));
        REQUIRE(top2.fetchParticles().at(1) == particles.at(2));
    }

    SECTION("Chain split up") {
        std::size_t n_chain_elements = 50;
        auto &toptypes = ctx.topologyRegistry();

        toptypes.addType("TA");

        ctx.boxSize() = {{10, 10, 10}};
        std::vector<readdy::model::Particle> topologyParticles;
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
            while (it2 != topology->graph().vertices().end()) {
                topology->addEdge(it.persistent_index(), it2.persistent_index());
                std::advance(it, 1);
                std::advance(it2, 1);
            }
        }

        {
            auto reactionFunction = [&](model::top::GraphTopology &top) {
                model::top::reactions::Recipe recipe(top);
                auto &vertices = top.graph().vertices();
                auto current_n_vertices = vertices.size();
                if (current_n_vertices > 1) {
                    auto edge = readdy::model::rnd::uniform_int<>(0, static_cast<int>(current_n_vertices - 2));
                    auto it1 = vertices.begin();
                    auto it2 = ++vertices.begin();
                    for (int i = 0; i < edge; ++i) {
                        ++it1;
                        ++it2;
                    }
                    recipe.removeEdge(it1.persistent_index(), it2.persistent_index());
                }
                return recipe;
            };
            auto rateFunction = [](const model::top::GraphTopology &top) {
                return top.nParticles() > 1 ? top.nParticles() / 50. : 0;
            };
            model::top::reactions::StructuralTopologyReaction reaction{"r", reactionFunction, rateFunction};
            reaction.create_child_topologies_after_reaction();
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
            for (time = 1; time < n_time_steps; ++time) {
                integrator->perform();
                topReactions->perform();
                forces->perform();
                kernel->evaluateObservables(time);
            }
        }

        auto topologies = kernel->stateModel().getTopologies();
        for (const auto topPtr : topologies) {
            REQUIRE(topPtr->nParticles() == 1);
            REQUIRE(topPtr->graph().vertices().size() == 1);
        }
    }

    SECTION("Chain split up with decay") {
        std::size_t n_chain_elements = 50;
        auto &toptypes = ctx.topologyRegistry();
        toptypes.addType("TA");

        ctx.boxSize() = {{10, 10, 10}};
        std::vector<readdy::model::Particle> topologyParticles;
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
            while (it2 != topology->graph().vertices().end()) {
                topology->addEdge(it.persistent_index(), it2.persistent_index());
                std::advance(it, 1);
                std::advance(it2, 1);
            }
        }

        {
            // split reaction
            auto reactionFunction = [&](model::top::GraphTopology &top) {
                model::top::reactions::Recipe recipe(top);
                auto &vertices = top.graph().vertices();
                auto current_n_vertices = vertices.size();
                if (current_n_vertices > 1) {
                    auto edge = readdy::model::rnd::uniform_int<>(0, static_cast<int>(current_n_vertices - 2));
                    auto it1 = vertices.begin();
                    auto it2 = ++vertices.begin();
                    for (int i = 0; i < edge; ++i) {
                        ++it1;
                        ++it2;
                    }
                    recipe.removeEdge(it1.persistent_index(), it2.persistent_index());
                }

                return recipe;
            };
            auto rateFunction = [](const model::top::GraphTopology &top) {
                return top.nParticles() > 1 ? top.nParticles() / 50. : 0;
            };
            model::top::reactions::StructuralTopologyReaction reaction{"r", reactionFunction, rateFunction};
            reaction.create_child_topologies_after_reaction();

            toptypes.addStructuralReaction("TA", reaction);
        }
        {
            // decay reaction
            auto reactionFunction = [&](model::top::GraphTopology &top) {
                model::top::reactions::Recipe recipe(top);
                if (top.graph().vertices().size() == 1) {
                    recipe.changeParticleType(top.graph().vertices().begin().persistent_index(),
                                              kernel->context().particleTypes().idOf("A"));
                } else {
                    throw std::logic_error("this reaction should only be executed when there is exactly "
                                           "one particle in the topology");
                }
                return recipe;
            };
            auto rateFunction = [](const model::top::GraphTopology &top) {
                return top.nParticles() > 1 ? 0 : 1;
            };
            model::top::reactions::StructuralTopologyReaction reaction{"r2", reactionFunction, rateFunction};
            reaction.create_child_topologies_after_reaction();
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
            for (time = 1; time < n_time_steps; ++time) {
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
            model::Context ctx {};

            ctx.particleTypes().addTopologyType("X1", 0);
            ctx.particleTypes().addTopologyType("X2", 0);
            ctx.particleTypes().addTopologyType("Y", 0);
            ctx.particleTypes().addTopologyType("Z", 0);
            ctx.topologyRegistry().addType("T");
            ctx.topologyRegistry().addType("T2");

            ctx.topologyRegistry().configureBondPotential("Y", "Z", {0., .1});
            SECTION("Scalar Constant Rate") {
                ctx.topologyRegistry().addSpatialReaction("connect: T(X1) + T(X2) -> T2(Y--Z)", 1e10, 1.);
            }
            SECTION("Constant Rate via function") {
                ctx.topologyRegistry().addSpatialReaction("connect: T(X1) + T(X2) -> T2(Y--Z)",
                                                          std::move(constant_rate_function),
                                                          1.);
            }

            Simulation sim(std::move(kernel), ctx);

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
            REQUIRE(top->nParticles() == 2);
            auto particles = top->fetchParticles();
            REQUIRE((particles[0].type() == sim.context().particleTypes().idOf("Y") ||
                     particles[0].type() == sim.context().particleTypes().idOf("Z")));
            if (particles[0].type() == sim.context().particleTypes().idOf("Y")) {
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
            SECTION("Constant Scalar Rate") {
                ctx.topologyRegistry().addSpatialReaction("connect: T(X) + T(X) -> T2(Y--Z) [self=true]", 1e10, 1.);
            }
            SECTION("Constant Rate via Function") {
                ctx.topologyRegistry().addSpatialReaction("connect: T(X) + T(X) -> T2(Y--Z) [self=true]",
                                                          std::move(constant_rate_function), 1.);
            }

            Simulation sim(kernel->name(), ctx);

            auto p1 = sim.createTopologyParticle("X", {0, 0, -.01});
            auto p2 = sim.createTopologyParticle("X", {0, 0, 0});
            auto p3 = sim.createTopologyParticle("X", {0, 0, .01});
            auto t = sim.addTopology("T", {p1, p2, p3});
            t->addEdge({0}, {1});
            t->addEdge({1}, {2});

            sim.run(10, 1e-3);

            auto topologies = sim.currentTopologies();

            REQUIRE(topologies.size() == 1);
            auto top = topologies.at(0);
            REQUIRE(top->type() == sim.context().topologyRegistry().idOf("T2"));
            REQUIRE(top->nParticles() == 3);
            auto particles = top->fetchParticles();
            for (auto ptype : {"X", "Y", "Z"}) {
                auto itX = std::find_if(particles.begin(), particles.end(), [&](const auto &p) {
                    return p.type() == sim.context().particleTypes().idOf(ptype);
                });
                INFO("at least one " + std::string(ptype) + " particle should be contained");
                REQUIRE(itX != particles.end());
                particles.erase(itX);
            }
        }

        SECTION("TTFusionNetwork") {
            model::Context ctx;

            ctx.boxSize() = {{30, 30, 30}};

            ctx.particleTypes().addTopologyType("head", 0);
            ctx.particleTypes().addTopologyType("core", 0);
            ctx.particleTypes().addTopologyType("tail", 0);
            ctx.particleTypes().addTopologyType("link", 0);

            ctx.topologyRegistry().addType("polymer");

            ctx.topologyRegistry().configureBondPotential("head", "core", {.0, .01});
            ctx.topologyRegistry().configureBondPotential("core", "tail", {.0, .01});
            ctx.topologyRegistry().configureBondPotential("core", "core", {.0, .01});
            ctx.topologyRegistry().configureBondPotential("core", "link", {.0, .01});
            ctx.topologyRegistry().configureBondPotential("head", "link", {.0, .01});
            ctx.topologyRegistry().configureBondPotential("tail", "link", {.0, .01});
            ctx.topologyRegistry().configureBondPotential("link", "link", {.0, .01});

            SECTION("min. number of edges too large for full connection") {
                SECTION("Constant Scalar Rate") {
                    ctx.topologyRegistry().addSpatialReaction(
                            "connect: polymer(core) + polymer(core) -> polymer(link--link) [self=true, distance>16]",
                            1e10, 1.1);
                }
                SECTION("Constant Rate via Function") {
                    ctx.topologyRegistry().addSpatialReaction(
                            "connect: polymer(core) + polymer(core) -> polymer(link--link) [self=true, distance>16]",
                            std::move(constant_rate_function), 1.1);
                }


                Simulation sim(kernel->name(), ctx);

                auto p11 = sim.createTopologyParticle("head", {-1, 0, 0});
                auto p12 = sim.createTopologyParticle("core", {0, 0, 0});
                auto p13 = sim.createTopologyParticle("core", {1, 0, 0});
                auto p14 = sim.createTopologyParticle("core", {2, 0, 0});
                auto p15 = sim.createTopologyParticle("core", {3, 0, 0});
                auto p16 = sim.createTopologyParticle("tail", {4, 0, 0});

                auto p21 = sim.createTopologyParticle("head", {3, 1, 1});
                auto p22 = sim.createTopologyParticle("core", {3, 0, 1});
                auto p23 = sim.createTopologyParticle("core", {3, -1, 1});
                auto p24 = sim.createTopologyParticle("core", {3, -2, 1});
                auto p25 = sim.createTopologyParticle("core", {3, -3, 1});
                auto p26 = sim.createTopologyParticle("tail", {3, -4, 1});

                auto p31 = sim.createTopologyParticle("head", {4, -3, 0});
                auto p32 = sim.createTopologyParticle("core", {3, -3, 0});
                auto p33 = sim.createTopologyParticle("core", {2, -3, 0});
                auto p34 = sim.createTopologyParticle("core", {1, -3, 0});
                auto p35 = sim.createTopologyParticle("core", {0, -3, 0});
                auto p36 = sim.createTopologyParticle("tail", {-1, -3, 0});

                auto p41 = sim.createTopologyParticle("head", {0, -4, 1});
                auto p42 = sim.createTopologyParticle("core", {0, -3, 1});
                auto p43 = sim.createTopologyParticle("core", {0, -2, 1});
                auto p44 = sim.createTopologyParticle("core", {0, -1, 1});
                auto p45 = sim.createTopologyParticle("core", {0, 0, 1});
                auto p46 = sim.createTopologyParticle("tail", {0, 1, 1});

                auto t1 = sim.addTopology("polymer", {p11, p12, p13, p14, p15, p16});
                auto t2 = sim.addTopology("polymer", {p21, p22, p23, p24, p25, p26});
                auto t3 = sim.addTopology("polymer", {p31, p32, p33, p34, p35, p36});
                auto t4 = sim.addTopology("polymer", {p41, p42, p43, p44, p45, p46});
                for (unsigned i = 0; i < 5; i++) {
                    t1->addEdge({i}, {i + 1});
                    t2->addEdge({i}, {i + 1});
                    t3->addEdge({i}, {i + 1});
                    t4->addEdge({i}, {i + 1});
                }
                sim.run(3, 1e-3);

                auto topologies = sim.currentTopologies();
                REQUIRE(topologies.size() == 1);

                auto top = topologies.at(0);
                auto particles = sim.currentParticles();
                unsigned n_links = 0;
                auto typeLink = ctx.particleTypes().infoOf("link").typeId;
                for (const auto & vref : top->graph()) {
                    auto pidx = vref.data().particleIndex;
                    if (particles.at(pidx).type() == typeLink) n_links++;
                }
                REQUIRE(6 == n_links);

                //for (auto p: particles) {
                //  if (p.type() == typeLink) n_links++;
                //}
                //REQUIRE( 6 == n_links );
            }
            SECTION("full connection") {
                SECTION("Constant Scalar Rate") {
                    ctx.topologyRegistry().addSpatialReaction(
                            "connect: polymer(core) + polymer(core) -> polymer(link--link) [self=true, distance>14]", // >14
                            1e10, 1.01);
                }
                SECTION("Constant Rate via Function") {
                    ctx.topologyRegistry().addSpatialReaction(
                            "connect: polymer(core) + polymer(core) -> polymer(link--link) [self=true, distance>14]", // >14
                            std::move(constant_rate_function), 1.01);
                }

                Simulation sim(kernel->name(), ctx);

                auto p11 = sim.createTopologyParticle("head", {-1, 0, 0});
                auto p12 = sim.createTopologyParticle("core", {0, 0, 0});
                auto p13 = sim.createTopologyParticle("core", {1, 0, 0});
                auto p14 = sim.createTopologyParticle("core", {2, 0, 0});
                auto p15 = sim.createTopologyParticle("core", {3, 0, 0});
                auto p16 = sim.createTopologyParticle("tail", {4, 0, 0});

                auto p21 = sim.createTopologyParticle("head", {3, 1, 1});
                auto p22 = sim.createTopologyParticle("core", {3, 0, 1});
                auto p23 = sim.createTopologyParticle("core", {3, -1, 1});
                auto p24 = sim.createTopologyParticle("core", {3, -2, 1});
                auto p25 = sim.createTopologyParticle("core", {3, -3, 1});
                auto p26 = sim.createTopologyParticle("tail", {3, -4, 1});

                auto p31 = sim.createTopologyParticle("head", {4, -3, 0});
                auto p32 = sim.createTopologyParticle("core", {3, -3, 0});
                auto p33 = sim.createTopologyParticle("core", {2, -3, 0});
                auto p34 = sim.createTopologyParticle("core", {1, -3, 0});
                auto p35 = sim.createTopologyParticle("core", {0, -3, 0});
                auto p36 = sim.createTopologyParticle("tail", {-1, -3, 0});

                auto p41 = sim.createTopologyParticle("head", {0, -4, 1});
                auto p42 = sim.createTopologyParticle("core", {0, -3, 1});
                auto p43 = sim.createTopologyParticle("core", {0, -2, 1});
                auto p44 = sim.createTopologyParticle("core", {0, -1, 1});
                auto p45 = sim.createTopologyParticle("core", {0, 0, 1});
                auto p46 = sim.createTopologyParticle("tail", {0, 1, 1});


                auto t1 = sim.addTopology("polymer", {p11, p12, p13, p14, p15, p16});
                auto t2 = sim.addTopology("polymer", {p21, p22, p23, p24, p25, p26});
                auto t3 = sim.addTopology("polymer", {p31, p32, p33, p34, p35, p36});
                auto t4 = sim.addTopology("polymer", {p41, p42, p43, p44, p45, p46});
                for (unsigned i = 0; i < 5; i++) {
                    t1->addEdge({i}, {i+1});
                    t2->addEdge({i}, {i+1});
                    t3->addEdge({i}, {i+1});
                    t4->addEdge({i}, {i+1});
                }
                sim.run(10, 1e-3);

                auto topologies = sim.currentTopologies();
                REQUIRE(topologies.size() == 1);

                auto top = topologies.at(0);
                auto particles = sim.currentParticles();
                unsigned n_links = 0;
                auto typeLink = ctx.particleTypes().infoOf("link").typeId;

                for (auto vref = top->graph().begin(); vref != top->graph().end(); vref++) {
                    auto pidx = vref->data().particleIndex;
                    if (vref->neighbors().empty()) continue;
                    if (particles.at(pidx).type() == typeLink) n_links++;
                }
                REQUIRE(8 == n_links);
            }
        }
    }
}
