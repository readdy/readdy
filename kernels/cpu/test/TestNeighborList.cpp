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

#include <gtest/gtest.h>
#include <readdy/api/SimulationLoop.h>
#include <readdy/kernel/cpu/CPUKernel.h>
#include <readdy/testing/NOOPPotential.h>


namespace cpu = readdy::kernel::cpu;
namespace m = readdy::model;

namespace {

using data_t = cpu::data::NLDataContainer;

struct TestNeighborList : ::testing::Test {

    std::unique_ptr<cpu::CPUKernel> kernel;
    readdy::ParticleTypeId typeIdA;

    std::unique_ptr<readdy::testing::NOOPPotentialOrder2> noop;

    TestNeighborList() : kernel(std::make_unique<cpu::CPUKernel>()) {
        auto &ctx = kernel->context();
        ctx.particleTypes().add("A", 1.);
        auto a_id = ctx.particleTypes()("A");
        readdy::scalar eductDistance = 1.2;
        ctx.reactions().addFusion("test", "A", "A", "A", 0., eductDistance);

        noop = std::make_unique<readdy::testing::NOOPPotentialOrder2>(a_id, a_id, 1.1, 0., 0.);
        ctx.potentials().addUserDefined(noop.get());
        typeIdA = ctx.particleTypes().idOf("A");
    }

};

class TestNeighborListImpl : public ::testing::TestWithParam<const char*> {};

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

TEST_F(TestNeighborList, ThreeBoxesNonPeriodic) {
    // maxcutoff is 1.2, system is 1.5 x 4 x 1.5, non-periodic, three cells
    auto &ctx = kernel->context();
    ctx.boxSize() = {{1.5, 4, 1.5}};
    ctx.periodicBoundaryConditions() = {{false, false, false}};

    readdy::kernel::cpu::thread_pool pool (readdy::readdy_default_n_threads());

    nl_t list(*kernel->getCPUKernelStateModel().getParticleData(), ctx, pool);

    auto &data = list.data();

    // Add three particles, two are in one outer box, the third on the other end and thus no neighbor
    const auto particles = std::vector<m::Particle>{
            m::Particle(0, -1.8, 0, typeIdA), m::Particle(0, -1.8, 0, typeIdA), m::Particle(0, 1.8, 0, typeIdA)
    };

    data.addParticles(particles);
    list.setUp(0, 1, {});

    EXPECT_TRUE(isPairInList(&list, 0, 1));
    EXPECT_TRUE(isPairInList(&list, 1, 0));
}

TEST_F(TestNeighborList, OneDirection) {
    // maxcutoff is 1.2, system is 4.8 x 5 x 5.1
    auto &ctx = kernel->context();
    ctx.boxSize() = {{1.2, 1.1, 2.8}};
    ctx.periodicBoundaryConditions() = {{false, false, true}};
    ctx.potentials().addBox("A", .0, {-.4, -.4, -1.3}, {.4, .4, 1.3});

    readdy::kernel::cpu::thread_pool pool (readdy::readdy_default_n_threads());
    nl_t list(*kernel->getCPUKernelStateModel().getParticleData(), ctx, pool);
    // Add three particles, one of which is in the neighborhood of the other two
    const auto particles = std::vector<m::Particle>{
            m::Particle(0, 0, -1.1, typeIdA), m::Particle(0, 0, .4, typeIdA), m::Particle(0, 0, 1.1, typeIdA)
    };
    std::vector<std::size_t> ids(particles.size());
    std::transform(particles.begin(), particles.end(), ids.begin(), [](const m::Particle& p) {return p.getId();});
    auto &data = list.data();
    data.addParticles(particles);

    list.setUp(0, 8, {});

    EXPECT_TRUE(isIdPairInList(&list, data, ids.at(0), ids.at(2)));
    EXPECT_TRUE(isIdPairInList(&list, data, ids.at(2), ids.at(0)));
    EXPECT_TRUE(isIdPairInList(&list, data, ids.at(1), ids.at(2)));
    EXPECT_TRUE(isIdPairInList(&list, data, ids.at(2), ids.at(1)));
    //EXPECT_FALSE(isIdPairInList(&list, data, ids.at(0), ids.at(1)));
    //EXPECT_FALSE(isIdPairInList(&list, data, ids.at(1), ids.at(0)));
}

TEST_F(TestNeighborList, AllNeighborsInCutoffSphere) {
    // maxcutoff is 1.2, system is 4 x 4 x 4, all directions periodic
    auto &ctx = kernel->context();
    ctx.boxSize() = {{4, 4, 4}};
    ctx.periodicBoundaryConditions() = {{true, true, true}};
    readdy::kernel::cpu::thread_pool pool (readdy::readdy_default_n_threads());
    nl_t list(*kernel->getCPUKernelStateModel().getParticleData(), ctx, pool);
    auto &data = list.data();
    // Create a few particles. In this box setup, all particles are neighbors.
    const auto particles = std::vector<m::Particle>{
            m::Particle(0, 0, 0, typeIdA), m::Particle(0, 0, 0, typeIdA), m::Particle(.3, 0, 0, typeIdA),
            m::Particle(0, .3, -.3, typeIdA), m::Particle(-.3, 0, .3, typeIdA), m::Particle(.3, -.3, 0, typeIdA)
    };

    data.addParticles(particles);
    list.setUp(0,1,{});
    for (size_t i = 0; i < 6; ++i) {
        for (size_t j = i + 1; j < 6; ++j) {
            EXPECT_TRUE(isPairInList(&list, i, j)) << "Particles " << i << " and " << j << " were not neighbors.";
            EXPECT_TRUE(isPairInList(&list, j, i)) << "Particles " << j << " and " << i << " were not neighbors.";
        }
    }
}


TEST(TestNeighborListImpl, DiffusionAndReaction) {
    using namespace readdy;
    std::unique_ptr<kernel::cpu::CPUKernel> kernel = std::make_unique<kernel::cpu::CPUKernel>();

    // A is absorbed and created by F, while the number of F stays constant, this test spans multiple timesteps
    kernel->context().particleTypes().add("A", 0.05);
    kernel->context().particleTypes().add("F", 0.0);
    kernel->context().particleTypes().add("V", 0.0);
    kernel->context().periodicBoundaryConditions() = {{true, true, true}};
    kernel->context().boxSize() = {{100, 10, 10}};

    const auto weightF = static_cast<readdy::scalar>(0.);
    const auto weightA = static_cast<readdy::scalar>(1.);
    kernel->context().reactions().addFusion("F+A->F", "A", "F", "F", .1, 2.0, weightF, weightA);
    //kernel->registerReaction<readdy::model::reactions::Fusion>("F+A->F2", "V", "V", "V", .1, 2.0, weightF, weightA);

    auto n3 = readdy::model::rnd::normal3<readdy::scalar>;
    // 120 F particles
    for (std::size_t i = 0; i < 100; ++i) {
        kernel->addParticle("F", n3(.0, 1.));
        kernel->addParticle("A", n3(0., 1.));
    }

    auto obs = kernel->observe().nParticles(1, std::vector<std::string>({"F", "A"}));
    obs->callback() = [&](const readdy::model::observables::NParticles::result_type &result) {
        EXPECT_EQ(result[0], 100);
    };
    auto connection = kernel->connectObservable(obs.get());

    {
        readdy::util::PerformanceNode pn("", false);
        readdy::api::SimulationLoop loop (kernel.get(), .01, pn);
        loop.useReactionScheduler("Gillespie");
        loop.skinSize() = .1;
        loop.run(100);
    }
}

TEST(TestNeighborListImpl, Diffusion) {
    using namespace readdy;
    std::unique_ptr<kernel::cpu::CPUKernel> kernel = std::make_unique<kernel::cpu::CPUKernel>();

    // A is absorbed and created by F, while the number of F stays constant, this test spans multiple timesteps
    auto& context = kernel->context();
    context.particleTypes().add("A", 0.05);
    context.particleTypes().add("F", 0.0);
    context.particleTypes().add("V", 0.0);

    scalar cutoff = 2.0;

    // just to have a cutoff
    context.reactions().addFusion("test", "V", "V", "V", .1, cutoff);
    context.periodicBoundaryConditions() = {{true, true, true}};
    context.boxSize() = {{100, 10, 10}};

    auto n3 = readdy::model::rnd::normal3<readdy::scalar>;
    // 120 F particles
    for (std::size_t i = 0; i < 50; ++i) {
        kernel->addParticle("F", n3(0., 1.));
        kernel->addParticle("A", n3(0., 1.));
    }
    auto obs = kernel->observe().nParticles(1);
    obs->callback() = (
            [&](const readdy::model::observables::NParticles::result_type &) {
                const auto neighbor_list = kernel->getCPUKernelStateModel().getNeighborList();

                for(std::size_t cell = 0; cell < neighbor_list->nCells(); ++cell) {
                    for(auto itParticle = neighbor_list->particlesBegin(cell);
                        itParticle != neighbor_list->particlesEnd(cell); ++itParticle) {
                        const auto &entry = neighbor_list->data().entry_at(*itParticle);
                        ASSERT_FALSE(entry.deactivated);

                        std::vector<std::size_t> neighbors;
                        neighbor_list->forEachNeighbor(*itParticle, [&](auto neighborIdx) {
                            const auto &neighborEntry = neighbor_list->data().entry_at(neighborIdx);
                            ASSERT_FALSE(neighborEntry.deactivated);
                            neighbors.push_back(neighborIdx);
                        });

                        std::size_t pidx = 0;
                        for(const auto &e : neighbor_list->data()) {
                            ASSERT_FALSE(e.deactivated);
                            if (pidx != *itParticle && bcs::dist(entry.pos, e.pos, context.boxSize(), context.periodicBoundaryConditions()) < cutoff) {
                                ASSERT_TRUE(std::find(neighbors.begin(), neighbors.end(), pidx) != neighbors.end());
                            }
                            ++pidx;
                        }
                    }
                }
            }
    );
    auto connection = kernel->connectObservable(obs.get());
    {
        readdy::util::PerformanceNode pn("", false);
        readdy::api::SimulationLoop loop(kernel.get(), .01, pn);
        loop.useReactionScheduler("Gillespie");
        loop.skinSize() = .1;
        loop.run(100);
    }
}

class TestCPUNeighborList : public ::testing::TestWithParam<std::array<readdy::scalar, 3>> {
public:

    auto getTestingBoxSize() const -> decltype(GetParam()) {
        return GetParam();
    }

};

TEST_P(TestCPUNeighborList, Diffuse) {
    using namespace readdy;

    using IndexPair = std::tuple<std::size_t, std::size_t>;

    kernel::cpu::CPUKernel kernel;

    auto& context = kernel.context();
    context.particleTypes().add("Test", 1.);
    auto id = context.particleTypes().idOf("Test");
    scalar cutoff = 4;
    context.reactions().addFusion("Fusion", id, id, id, .001, cutoff);
    context.boxSize() = getTestingBoxSize();
    /*context.boxSize()[0] = 10;
    context.boxSize()[1] = 5;
    context.boxSize()[2] = 5;*/
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
        kernel.stateModel().addParticle(particle);
    }

    kernel.initialize();

    kernel.stateModel().initializeNeighborList(0);

    auto integrator = kernel.actions().eulerBDIntegrator(.1);
    auto reactionHandler = kernel.actions().uncontrolledApproximation(.1);

    const auto &data = *kernel.getCPUKernelStateModel().getParticleData();

    // collect all pairs of particles that are closer than cutoff, these should (uniquely) be in the NL
    std::unordered_set<IndexPair, util::ForwardTupleHasher<IndexPair>, util::ForwardTupleEquality<IndexPair>> pairs;
    for(auto t = 0U; t < n_steps; ++t) {
        integrator->perform();

        kernel.stateModel().updateNeighborList();
        {
            pairs.clear();
            std::size_t ix1 = 0;
            for(const auto &e1 : data) {
                std::size_t ix2 = 0;
                for(const auto &e2 : data) {
                    if(ix1 != ix2 && !e1.deactivated && !e2.deactivated && bcs::dist(e1.pos, e2.pos, context.boxSize(), context.periodicBoundaryConditions()) < cutoff) {
                        pairs.insert(std::make_tuple(ix1, ix2));
                    }
                    ++ix2;
                }
                ++ix1;
            }
        }
        auto pairsCopy = pairs;

        const auto &neighborList = *kernel.getCPUKernelStateModel().getNeighborList();
        for (auto cell = 0U; cell < neighborList.nCells(); ++cell) {
            for(auto it = neighborList.particlesBegin(cell); it != neighborList.particlesEnd(cell); ++it) {
                const auto &entry = data.entry_at(*it);
                EXPECT_FALSE(entry.deactivated) << "A deactivated entry should not end up in the NL";
                neighborList.forEachNeighbor(*it, cell, [&](const std::size_t neighborIndex) {
                    const auto &neighbor = data.entry_at(neighborIndex);
                    EXPECT_FALSE(neighbor.deactivated) << "A deactivated entry should not end up in the NL";
                    // we got a pair
                    if(bcs::dist(entry.pos, neighbor.pos, context.boxSize(), context.periodicBoundaryConditions()) < cutoff) {
                        auto findIt = pairs.find(std::make_tuple(*it, neighborIndex));
                        std::stringstream ss;
                        if(pairsCopy.find(std::make_tuple(*it, neighborIndex)) != pairsCopy.end()) {
                            ss << "A particle pairing was more than once in the NL";
                        } else {
                            ss << "A particle pairing was not contained in the NL";
                        }

                        ASSERT_NE(findIt, pairs.end()) << ss.str();
                        if(findIt != pairs.end()) {
                            pairs.erase(findIt);
                        }
                    }
                });
            }
        }
        EXPECT_EQ(pairs.size(), 0) << "Some pairs were not contained in the NL";

        reactionHandler->perform();
    }
}

}

INSTANTIATE_TEST_CASE_P(
        CPUNeighborList,
        TestCPUNeighborList,
        ::testing::Values(
                std::array<readdy::scalar, 3>{{15, 15, 15}},
                std::array<readdy::scalar, 3>{{10, 10, 10}},
                std::array<readdy::scalar, 3>{{10, 5, 5}},
                std::array<readdy::scalar, 3>{{5, 5, 10}},
                std::array<readdy::scalar, 3>{{15, 5, 10}},
                std::array<readdy::scalar, 3>{{5, 10, 5}},
                std::array<readdy::scalar, 3>{{5, 5, 5}}
        )
);

