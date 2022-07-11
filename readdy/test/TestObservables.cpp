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
 * @file TestObservables.cpp
 * @brief Testing observables
 * @author clonker
 * @date 02.05.16
 */

#include <catch2/catch_template_test_macros.hpp>

#include <readdy/plugin/KernelProvider.h>
#include <readdy/api/Simulation.h>
#include <readdy/testing/KernelTest.h>
#include <readdy/testing/Utils.h>

namespace m = readdy::model;

using namespace readdytesting::kernel;

TEMPLATE_TEST_CASE("Test observables", "[observables]", SingleCPU, CPU) {
    auto kernel = create<TestType>();
    auto &context = kernel->context();
    auto &stateModel = kernel->stateModel();

    SECTION("Particle positions") {
        const unsigned int n_particles = 100;
        context.particleTypes().add("type", 1.);
        const readdy::scalar timeStep = 1.0;
        const auto particleTypeId = context.particleTypes().idOf("type");
        const auto particles = std::vector<m::Particle>(n_particles, m::Particle(0, 0, 0, particleTypeId));
        stateModel.addParticles(particles);
        auto &&obs = kernel->observe().positions(3);
        auto &&connection = kernel->connectObservable(obs.get());

        kernel->initialize();

        auto &&integrator = kernel->actions().createIntegrator("EulerBDIntegrator", timeStep);
        for (readdy::TimeStep t = 0; t < 100; t++) {
            integrator->perform();
            kernel->evaluateObservables(t);
        }

        const auto &result = obs->getResult();
        const auto &&positions = kernel->stateModel().getParticlePositions();
        auto it_pos = positions.begin();
        int j = 0;
        for (auto it = result.begin(); it != result.end(); it = std::next(it)) {
            REQUIRE(*it == *it_pos);
            it_pos++;
            ++j;
        }
        REQUIRE(j == 100);
        connection.disconnect();
    }

    SECTION("Topologies") {
        context.particleTypes().add("Topology A", 1.0, readdy::model::particleflavor::TOPOLOGY);
        context.particleTypes().add("Topology B", 1.0, readdy::model::particleflavor::TOPOLOGY);
        context.particleTypes().add("Topology Invalid Type", 1.0, readdy::model::particleflavor::TOPOLOGY);
        context.particleTypes().add("A", 1.0, readdy::model::particleflavor::NORMAL);

        context.topologyRegistry().configureBondPotential("Topology A", "Topology A", {10, 10});
        context.topologyRegistry().configureBondPotential("Topology A", "Topology B", {10, 10});
        context.topologyRegistry().configureBondPotential("Topology B", "Topology B", {10, 10});

        context.boxSize() = {{10, 10, 10}};

        std::size_t n_chain_elements = 50;
        auto &toptypes = context.topologyRegistry();
        toptypes.addType("TA");

        std::vector<readdy::model::Particle> topologyParticles;
        {
            topologyParticles.reserve(n_chain_elements);
            for (std::size_t i = 0; i < n_chain_elements; ++i) {
                const auto id = context.particleTypes().idOf("Topology A");
                topologyParticles.emplace_back(-5 + i * 10. / static_cast<readdy::scalar>(n_chain_elements), 0, 0, id);
            }
        }
        auto topology = kernel->stateModel().addTopology(toptypes.idOf("TA"), topologyParticles);
        {
            std::size_t it = 0;
            std::size_t it2 = 1;
            while (it2 < topology->graph().vertices().size()) {
                topology->addEdge({it}, {it2});
                ++it;
                ++it2;
            }
        }

        {
            // split reaction
            auto reactionFunction = [&](readdy::model::top::GraphTopology &top) {
                readdy::model::top::reactions::Recipe recipe(top);
                auto &vertices = top.graph().vertices();
                auto current_n_vertices = vertices.size();
                if (current_n_vertices > 1) {
                    auto edge = readdy::model::rnd::uniform_int<>(0, static_cast<int>(current_n_vertices - 2));
                    std::size_t it1 = 0;
                    std::size_t it2 = 1;
                    for (int i = 0; i < edge; ++i) {
                        ++it1;
                        ++it2;
                    }
                    recipe.removeEdge({it1}, {it2});
                }

                return recipe;
            };
            auto rateFunction = [](const readdy::model::top::GraphTopology &top) {
                return top.graph().nVertices() > 1 ? top.graph().nVertices() / 50. : 0;
            };
            readdy::model::top::reactions::StructuralTopologyReaction reaction{"r", reactionFunction, rateFunction};
            reaction.create_child_topologies_after_reaction();

            toptypes.addStructuralReaction("TA", reaction);
        }
        {
            // decay reaction
            auto reactionFunction = [&](readdy::model::top::GraphTopology &top) {
                readdy::model::top::reactions::Recipe recipe(top);
                if (top.graph().vertices().size() == 1) {
                    recipe.changeParticleType(top.graph().begin().persistent_index(), context.particleTypes().idOf("A"));
                } else {
                    throw std::logic_error("this reaction should only be executed when there is exactly "
                                           "one particle in the topology");
                }
                return recipe;
            };
            auto rateFunction = [](const readdy::model::top::GraphTopology &top) {
                return top.nParticles() > 1 ? 0 : 1;
            };
            readdy::model::top::reactions::StructuralTopologyReaction reaction{"r2", reactionFunction, rateFunction};
            reaction.create_child_topologies_after_reaction();
            toptypes.addStructuralReaction("TA", reaction);
        }

        {
            auto integrator = kernel->actions().createIntegrator("EulerBDIntegrator", 1.0);
            auto forces = kernel->actions().calculateForces();
            auto topReactions = kernel->actions().evaluateTopologyReactions(1.0);

            std::size_t time = 0;
            std::size_t n_time_steps = 500;

            auto obs = kernel->observe().topologies(1);
            obs->setCallback([&](const readdy::model::observables::Topologies::result_type &value) {
                auto tops = kernel->stateModel().getTopologies();
                REQUIRE(value.size() == tops.size());
                for (auto its = std::make_pair(tops.begin(), value.begin());
                     its.first != tops.end(); ++its.first, ++its.second) {
                    auto topPtr = *its.first;
                    const auto &record = *its.second;
                    const auto &topParticles = topPtr->particleIndices();
                    const auto &recordParticles = record.particleIndices;
                    auto contains1 = std::all_of(recordParticles.begin(), recordParticles.end(), [&](auto idx) {
                        return std::find(topParticles.begin(), topParticles.end(), idx) != topParticles.end();
                    });
                    auto contains2 = std::all_of(topParticles.begin(), topParticles.end(), [&](auto idx) {
                        return std::find(recordParticles.begin(), recordParticles.end(), idx) != recordParticles.end();
                    });
                    // record.particleIndices should be contained in topology particles
                    REQUIRE(contains1);
                    // topology particles should be contained in record.particleIndices
                    REQUIRE(contains2);

                    //auto topEdges = topPtr->graph().edges();
                    std::vector<readdy::model::top::Graph::Edge> topEdges;
                    topPtr->graph().findEdges([&topEdges](const auto &e) { topEdges.push_back(e); });
                    REQUIRE(topEdges.size() == record.edges.size());

                    contains1 = std::all_of(record.edges.begin(), record.edges.end(), [&](const auto &edge) {
                        std::size_t ix1 = std::get<0>(edge);
                        std::size_t ix2 = std::get<1>(edge);
                        for (const auto &topEdge : topEdges) {
                            const auto &[v1, v2] = topEdge;
                            if (v1.value == ix1 && v2.value == ix2) {
                                return true;
                            }
                            if (v1.value == ix2 && v2.value == ix1) {
                                return true;
                            }
                        }
                        return false;
                    });

                    contains2 = std::all_of(topEdges.begin(), topEdges.end(), [&](const auto &e) {
                        auto vtup1 = std::make_tuple(std::get<0>(e).value, std::get<1>(e).value);
                        auto vtup2 = std::make_tuple(std::get<1>(e).value, std::get<0>(e).value);

                        auto find1 = std::find(record.edges.begin(), record.edges.end(), vtup1);
                        auto find2 = std::find(record.edges.begin(), record.edges.end(), vtup2);

                        return find1 != record.edges.end() || find2 != record.edges.end();
                    });

                    // there are no non existent edges in the record
                    REQUIRE(contains1);
                    // all edges are in the record
                    REQUIRE(contains2);
                }
            });
            auto connection = kernel->connectObservable(obs.get());

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
    }
    SECTION("Forces") {
        // Setup particles
        context.particleTypes().add("A", 42.);
        context.particleTypes().add("B", 1337.);
        const auto typeIdA = context.particleTypes().idOf("A");
        const auto typeIdB = context.particleTypes().idOf("B");
        const unsigned int n_particles = 2; // There will be 55 Bs
        const auto particlesA = std::vector<m::Particle>(n_particles, m::Particle(0, 0, 0, typeIdA));
        const auto particlesB = std::vector<m::Particle>(n_particles + 5, m::Particle(0, 0, 0, typeIdB));
        kernel->stateModel().addParticles(particlesA);
        kernel->stateModel().addParticles(particlesB);
        {
            // Check if result has correct size
            // Check that empty particleType argument gives correct object, namely all forces
            auto &&obsA = kernel->observe().forces(1, std::vector<std::string>{"A"});
            auto &&obsB = kernel->observe().forces(1, std::vector<std::string>{"B"});
            auto &&obsBoth = kernel->observe().forces(1);
            auto &&connectionA = kernel->connectObservable(obsA.get());
            auto &&connectionB = kernel->connectObservable(obsB.get());
            auto &&connectionBoth = kernel->connectObservable(obsBoth.get());
            // Evaluate twice to ensure that results do not accumulate
            kernel->evaluateObservables(0);
            kernel->evaluateObservables(1);
            const auto &resA = obsA->getResult();
            const auto &resB = obsB->getResult();
            const auto &resBoth = obsBoth->getResult();
            REQUIRE(resA.size() == n_particles);
            REQUIRE(resB.size() == n_particles+5);
            REQUIRE(resBoth.size() == n_particles + n_particles + 5);
            readdy::Vec3 zero = readdy::Vec3(0, 0, 0);
            for (auto force : resBoth) {
                REQUIRE(force == zero);
            }
        }
        // Two particles C and C with radius 1 and harmonic repulsion at distance 1.5 -> force = kappa * (radiiSum - 1.5)
        context.periodicBoundaryConditions() = {{false, false, false}};
        context.boxSize() = {{5, 5, 5}};
        context.particleTypes().add("C", 1.);
        context.potentials().addBox("A", .001, {-2.4, -2.4, -2.4}, {4.8, 4.8, 4.8});
        context.potentials().addBox("B", .001, {-2.4, -2.4, -2.4}, {4.8, 4.8, 4.8});
        context.potentials().addBox("C", .001, {-2.4, -2.4, -2.4}, {4.8, 4.8, 4.8});
        const auto typeIdC = context.particleTypes().idOf("C");
        const auto particlesC = std::vector<m::Particle>{m::Particle(0, 0, 0, typeIdC), m::Particle(0, -1.5, 0, typeIdC)};
        kernel->stateModel().addParticles(particlesC);

        context.potentials().addHarmonicRepulsion("C", "C", 2.0, 2.0);

        auto &&initNeighborList = kernel->actions().createNeighborList(context.calculateMaxCutoff());
        auto &&nl = kernel->actions().updateNeighborList();
        auto &&forces = kernel->actions().calculateForces();
        kernel->initialize();
        {
            auto obsC = kernel->observe().forces(1, std::vector<std::string>{"C"});
            auto connectionC = kernel->connectObservable(obsC.get());
            initNeighborList->perform();
            nl->perform();
            forces->perform();
            kernel->evaluateObservables(2);
            const auto &resC = obsC->getResult();
            readdy::Vec3 force0 = readdy::Vec3(0., 1., 0.);
            readdy::Vec3 force1 = readdy::Vec3(0., -1., 0.);
            REQUIRE(resC.size() == 2);
            REQUIRE((resC[0] == force0 || resC[1] == force0));
            REQUIRE((resC[1] == force1 || resC[0] == force1));
        }
    }
}
