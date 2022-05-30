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
 * Test different contexts w.r.t. boxsize and periodicity, perform setupBoxes() and see if that worked.
 * Then add a few particles and perform fillBoxes(). The result of this is the actual neighborlist, which is
 * checked as well.
 *
 * @file CPUTestNeighborList.cpp
 * @brief Test the neighborlist object of the CPU kernel.
 * @author chrisfroe
 * @author clonker
 * @date 23.08.16
 */

#include <cmath>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <readdy/api/SimulationLoop.h>
#include <readdy/kernel/cpu/CPUKernel.h>
#include <readdy/testing/NOOPPotential.h>


namespace cpu = readdy::kernel::cpu;
namespace m = readdy::model;

using data_t = cpu::data::NLDataContainer;
using nl_t = readdy::kernel::cpu::nl::CompactCellLinkedList;

auto isPairInList = [](nl_t *pairs, std::size_t idx1, std::size_t idx2) {
    bool foundOneDirection {false};
    bool foundOtherDirection {false};
    {
        pairs->forEachNeighbor(idx1, [&](auto neighborIdx) {
            if(neighborIdx == idx2) {
                foundOneDirection = true;
            }
        });
    }
    {
        pairs->forEachNeighbor(idx2, [&](auto neighborIdx) {
            if(neighborIdx == idx1) {
                foundOtherDirection = true;
            }
        });
    }
    return foundOneDirection && foundOtherDirection;
};

auto isIdPairInList = [](nl_t *pairs, readdy::kernel::cpu::data::DefaultDataContainer &data, std::size_t id1, std::size_t id2) {
    return isPairInList(pairs, data.getIndexForId(id1), data.getIndexForId(id2));
};


TEST_CASE("Test cpu neighbor list", "[cpu]") {
    using namespace readdy;

    auto kernel = std::make_unique<cpu::CPUKernel>();
    auto &ctx = kernel->context();

    SECTION("Basic behavior") {
        ctx.particleTypes().add("A", 1.);
        auto a_id = ctx.particleTypes()("A");
        readdy::scalar eductDistance = 1.2;
        ctx.reactions().addFusion("test", "A", "A", "A", 0., eductDistance);

        auto noop = std::make_unique<readdy::testing::NOOPPotentialOrder2>(a_id, a_id, 1.1, 0., 0.);
        ctx.potentials().addUserDefined(noop.get());
        auto typeIdA = ctx.particleTypes().idOf("A");

        SECTION("Three boxes and non periodic") {
            ctx.boxSize() = {{1.5, 4, 1.5}};
            ctx.periodicBoundaryConditions() = {{false, false, false}};

            readdy::kernel::cpu::thread_pool pool(readdy::readdy_default_n_threads());

            nl_t list(*kernel->getCPUKernelStateModel().getParticleData(), ctx, pool);

            auto &data = list.data();

            // Add three particles, two are in one outer box, the third on the other end and thus no neighbor
            const auto particles = std::vector<m::Particle>{
                    m::Particle(0, -1.8, 0, typeIdA), m::Particle(0, -1.8, 0, typeIdA), m::Particle(0, 1.8, 0, typeIdA)
            };

            data.addParticles(particles);
            list.setUp(ctx.calculateMaxCutoff(), 1);

            REQUIRE(isPairInList(&list, 0, 1));
            REQUIRE(isPairInList(&list, 1, 0));
        }

        SECTION("Periodic in one direction") {
            // maxcutoff is 1.2, system is 4.8 x 5 x 5.1
            ctx.boxSize() = {{1.2, 1.1, 2.8}};
            ctx.periodicBoundaryConditions() = {{false, false, true}};
            ctx.potentials().addBox("A", .0, {-.4, -.4, -1.3}, {.4, .4, 1.3});

            readdy::kernel::cpu::thread_pool pool(readdy::readdy_default_n_threads());
            nl_t list(*kernel->getCPUKernelStateModel().getParticleData(), ctx, pool);
            // Add three particles, one of which is in the neighborhood of the other two
            const auto particles = std::vector<m::Particle>{
                    m::Particle(0, 0, -1.1, typeIdA), m::Particle(0, 0, .4, typeIdA), m::Particle(0, 0, 1.1, typeIdA)
            };
            std::vector<std::size_t> ids(particles.size());
            std::transform(particles.begin(), particles.end(), ids.begin(),
                           [](const m::Particle &p) { return p.id(); });
            auto &data = list.data();
            data.addParticles(particles);

            list.setUp(ctx.calculateMaxCutoff(), 8);

            REQUIRE(isIdPairInList(&list, data, ids.at(0), ids.at(2)));
            REQUIRE(isIdPairInList(&list, data, ids.at(2), ids.at(0)));
            REQUIRE(isIdPairInList(&list, data, ids.at(1), ids.at(2)));
            REQUIRE(isIdPairInList(&list, data, ids.at(2), ids.at(1)));
        }
        SECTION("All neighbors are within cutoff, periodic") {
            // maxcutoff is 1.2, system is 4 x 4 x 4, all directions periodic
            auto &ctx = kernel->context();
            ctx.boxSize() = {{4, 4, 4}};
            ctx.periodicBoundaryConditions() = {{true, true, true}};
            readdy::kernel::cpu::thread_pool pool(readdy::readdy_default_n_threads());
            nl_t list(*kernel->getCPUKernelStateModel().getParticleData(), ctx, pool);
            auto &data = list.data();
            // Create a few particles. In this box setup, all particles are neighbors.
            const auto particles = std::vector<m::Particle>{
                    m::Particle(0, 0, 0, typeIdA), m::Particle(0, 0, 0, typeIdA), m::Particle(.3, 0, 0, typeIdA),
                    m::Particle(0, .3, -.3, typeIdA), m::Particle(-.3, 0, .3, typeIdA), m::Particle(.3, -.3, 0, typeIdA)
            };

            data.addParticles(particles);
            list.setUp(ctx.calculateMaxCutoff(), 1);
            for (size_t i = 0; i < 6; ++i) {
                for (size_t j = i + 1; j < 6; ++j) {
                    REQUIRE(isPairInList(&list, i, j));
                    REQUIRE(isPairInList(&list, j, i));
                }
            }
        }
    }

    SECTION("Diffusion and reaction") {
        // A is absorbed and created by F, while the number of F stays constant, this test spans multiple timesteps
        ctx.particleTypes().add("A", 0.05);
        ctx.particleTypes().add("F", 0.0);
        ctx.particleTypes().add("V", 0.0);
        ctx.periodicBoundaryConditions() = {{true, true, true}};
        ctx.boxSize() = {{100, 10, 10}};

        const auto weightF = static_cast<readdy::scalar>(0.);
        const auto weightA = static_cast<readdy::scalar>(1.);
        ctx.reactions().addFusion("F+A->F", "A", "F", "F", .1, 2.0, weightF, weightA);

        auto n3 = readdy::model::rnd::normal3<readdy::scalar>;
        // 120 F particles
        for (std::size_t i = 0; i < 100; ++i) {
            kernel->addParticle("F", n3(.0, 1.));
            kernel->addParticle("A", n3(0., 1.));
        }

        auto obs = kernel->observe().nParticles(1, std::vector<std::string>({"F", "A"}));
        obs->setCallback([&](const readdy::model::observables::NParticles::result_type &result) {
            REQUIRE(result[0] == 100);
        });
        auto connection = kernel->connectObservable(obs.get());

        {
            readdy::api::SimulationLoop loop(kernel.get(), .01);
            loop.useReactionScheduler("Gillespie");
            loop.neighborListCutoff() += 0.1;
            loop.run(100);
        }
    }
    SECTION("Diffusion") {
        // A is absorbed and created by F, while the number of F stays constant, this test spans multiple timesteps
        ctx.particleTypes().add("A", 0.05);
        ctx.particleTypes().add("F", 0.0);
        ctx.particleTypes().add("V", 0.0);

        auto cutoff = 2.0;

        // just to have a cutoff
        ctx.reactions().addFusion("test", "V", "V", "V", .1, cutoff);
        ctx.periodicBoundaryConditions() = {{true, true, true}};
        ctx.boxSize() = {{100, 10, 10}};

        auto n3 = readdy::model::rnd::normal3<readdy::scalar>;
        // 120 F particles
        for (std::size_t i = 0; i < 50; ++i) {
            kernel->addParticle("F", n3(0., 1.));
            kernel->addParticle("A", n3(0., 1.));
        }
        auto obs = kernel->observe().nParticles(1);
        obs->setCallback(
                [&](const readdy::model::observables::NParticles::result_type &) {
                    const auto neighbor_list = kernel->getCPUKernelStateModel().getNeighborList();

                    for(std::size_t cell = 0; cell < neighbor_list->nCells(); ++cell) {
                        for(auto itParticle = neighbor_list->particlesBegin(cell);
                            itParticle != neighbor_list->particlesEnd(cell); ++itParticle) {
                            const auto &entry = neighbor_list->data().entry_at(*itParticle);
                            REQUIRE_FALSE(entry.deactivated);

                            std::vector<std::size_t> neighbors;
                            neighbor_list->forEachNeighbor(*itParticle, [&](auto neighborIdx) {
                                const auto &neighborEntry = neighbor_list->data().entry_at(neighborIdx);
                                REQUIRE_FALSE(neighborEntry.deactivated);
                                neighbors.push_back(neighborIdx);
                            });

                            std::size_t pidx = 0;
                            for(const auto &e : neighbor_list->data()) {
                                REQUIRE_FALSE(e.deactivated);
                                if (pidx != *itParticle && bcs::dist(entry.pos, e.pos, ctx.boxSize(),
                                        ctx.periodicBoundaryConditions()) < cutoff) {
                                    REQUIRE(std::find(neighbors.begin(), neighbors.end(), pidx) != neighbors.end());
                                }
                                ++pidx;
                            }
                        }
                    }
                }
        );
        auto connection = kernel->connectObservable(obs.get());
        {
            readdy::api::SimulationLoop loop(kernel.get(), .01);
            loop.useReactionScheduler("Gillespie");
            loop.neighborListCutoff() += 0.1;
            loop.run(100);
        }
    }
    SECTION("Diffusion in different boxes") {
        auto boxSize = GENERATE(std::array<readdy::scalar, 3>{{15, 15, 15}},
                                std::array<readdy::scalar, 3>{{10, 10, 10}},
                                std::array<readdy::scalar, 3>{{10, 5, 5}},
                                std::array<readdy::scalar, 3>{{5, 5, 10}},
                                std::array<readdy::scalar, 3>{{15, 5, 10}},
                                std::array<readdy::scalar, 3>{{5, 10, 5}},
                                std::array<readdy::scalar, 3>{{5, 5, 5}});

        INFO(fmt::format("Testing diffusion in box size ({}, {}, {})", boxSize[0], boxSize[1], boxSize[2]));

        using IndexPair = std::tuple<std::size_t, std::size_t>;

        auto& context = kernel->context();
        context.particleTypes().add("Test", 1.);
        auto id = context.particleTypes().idOf("Test");
        scalar cutoff = 4;
        context.reactions().addFusion("Fusion", id, id, id, .001, cutoff);
        context.boxSize() = boxSize;
        bool periodic = true;
        context.periodicBoundaryConditions()[0] = periodic;
        context.periodicBoundaryConditions()[1] = periodic;
        context.periodicBoundaryConditions()[2] = periodic;

#ifdef READDY_DEBUG
        auto n_steps = 3U;
        auto n_particles = 500;
#else
        auto n_steps = 6U;
    auto n_particles = 1000;
#endif

        for(auto i = 0; i < n_particles; ++i) {
            model::Particle particle(model::rnd::uniform_real<scalar>(-.5*context.boxSize()[0], .5*context.boxSize()[0]),
                                     model::rnd::uniform_real<scalar>(-.5*context.boxSize()[1], .5*context.boxSize()[1]),
                                     model::rnd::uniform_real<scalar>(-.5*context.boxSize()[2], .5*context.boxSize()[2]), id);
            kernel->stateModel().addParticle(particle);
        }

        kernel->initialize();

        kernel->stateModel().initializeNeighborList(context.calculateMaxCutoff());

        auto integrator = kernel->actions().eulerBDIntegrator(.1);
        auto reactionHandler = kernel->actions().uncontrolledApproximation(.1);

        const auto &data = *kernel->getCPUKernelStateModel().getParticleData();

        // collect all pairs of particles that are closer than cutoff, these should (uniquely) be in the NL
        std::unordered_set<IndexPair, util::ForwardTupleHasher<IndexPair>, util::ForwardTupleEquality<IndexPair>> pairs;
        for(auto t = 0U; t < n_steps; ++t) {
            integrator->perform();

            kernel->stateModel().updateNeighborList();
            {
                pairs.clear();
                std::size_t ix1 = 0;
                for(const auto &e1 : data) {
                    std::size_t ix2 = 0;
                    for(const auto &e2 : data) {
                        if(ix1 != ix2 && !e1.deactivated && !e2.deactivated && bcs::dist(e1.pos, e2.pos,
                                context.boxSize(), context.periodicBoundaryConditions()) < cutoff) {
                            pairs.insert(std::make_tuple(ix1, ix2));
                        }
                        ++ix2;
                    }
                    ++ix1;
                }
            }
            auto pairsCopy = pairs;

            const auto &neighborList = *kernel->getCPUKernelStateModel().getNeighborList();
            for (auto cell = 0U; cell < neighborList.nCells(); ++cell) {
                for(auto it = neighborList.particlesBegin(cell); it != neighborList.particlesEnd(cell); ++it) {
                    const auto &entry = data.entry_at(*it);
                    REQUIRE_FALSE(entry.deactivated);
                    neighborList.forEachNeighbor(*it, cell, [&](const std::size_t neighborIndex) {
                        const auto &neighbor = data.entry_at(neighborIndex);
                        REQUIRE_FALSE(neighbor.deactivated);
                        // we got a pair
                        if(bcs::dist(entry.pos, neighbor.pos, context.boxSize(),
                                context.periodicBoundaryConditions()) < cutoff) {
                            auto findIt = pairs.find(std::make_tuple(*it, neighborIndex));

                            REQUIRE(findIt != pairs.end());
                            if(findIt != pairs.end()) {
                                pairs.erase(findIt);
                            }
                        }
                    });
                }
            }
            REQUIRE(pairs.empty());

            reactionHandler->perform();
        }
    }
}
