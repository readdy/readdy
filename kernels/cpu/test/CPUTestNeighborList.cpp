/**
 * Test different contexts w.r.t. boxsize and periodicity, perform setupBoxes() and see if that worked.
 * Then add a few particles and perform fillBoxes(). The result of this is the actual neighborlist, which is
 * checked as well.
 *
 * @file CPUTestNeighborList.cpp
 * @brief Test the neighborlist object of the CPU kernel.
 * @author chrisfroe
 * @date 23.08.16
 */

#include <gtest/gtest.h>
#include <readdy/kernel/cpu/model/NeighborList.h>
#include <readdy/kernel/cpu/CPUKernel.h>
#include <readdy/kernel/singlecpu/SingleCPUKernel.h>
#include <readdy/testing/NOOPPotential.h>

namespace cpu = readdy::kernel::cpu;
namespace scpu = readdy::kernel::singlecpu;
namespace cpum = cpu::model;
namespace scpum = scpu::model;
namespace m = readdy::model;

namespace {

struct TestNeighborList : ::testing::Test {

    std::unique_ptr<cpu::CPUKernel> kernel;
    unsigned int typeIdA;

    TestNeighborList() : kernel(std::make_unique<cpu::CPUKernel>()) {
        auto &ctx = kernel->getKernelContext();
        ctx.setDiffusionConstant("A", 1.0);
        double eductDistance = 1.2;
        ctx.registerReaction(kernel->createFusionReaction("test", "A", "A", "A", 0., eductDistance));

        ctx.registerPotential(std::make_unique<readdy::testing::NOOPPotentialOrder2>(1.1, 0., 0.), "A", "A");
        typeIdA = ctx.getParticleTypeID("A");
        ctx.configure();
    }

};

auto isPairInList = [](readdy::kernel::cpu::model::NeighborList::container_t *pairs, unsigned long idx1,
                       unsigned long idx2) {
    using neighbor_t = readdy::kernel::cpu::model::NeighborList::neighbor_t;
    const auto neighborsIt = pairs->find(idx1);
    if (neighborsIt != pairs->end()) {
        const auto neighbors = neighborsIt->second;
        return std::find_if(neighbors.begin(), neighbors.end(), [idx2](const neighbor_t &neighbor) {
            return neighbor.idx == idx2;
        }) != neighbors.end();
    }
    return false;
};

auto getNumberPairs = [](readdy::kernel::cpu::model::NeighborList::container_t *pairs) {
    auto pairsIt = pairs->cbegin();
    size_t numberPairs = 0;
    while (pairsIt != pairs->cend()) {
        numberPairs += pairsIt->second.size();
        ++pairsIt;
    }
    return numberPairs;
};

TEST_F(TestNeighborList, ThreeBoxesNonPeriodic) {
    // maxcutoff is 1.2, system is 1.5 x 4 x 1.5, non-periodic, three boxes
    auto &ctx = kernel->getKernelContext();
    ctx.setBoxSize(1.5, 4, 1.5);
    ctx.setPeriodicBoundary(false, false, false);
    readdy::kernel::cpu::util::Config conf;
    cpum::NeighborList list(&ctx, &conf);
    list.setupBoxes();
    // Add three particles, two are in one outer box, the third on the other end and thus no neighbor
    const auto particles = std::vector<m::Particle>{
            m::Particle(0, -1.8, 0, typeIdA), m::Particle(0, -1.8, 0, typeIdA), m::Particle(0, 1.8, 0, typeIdA)
    };
    scpum::ParticleData data;
    data.addParticles(particles);
    list.fillBoxes(data);
    auto pairs = list.pairs.get();
    EXPECT_EQ(getNumberPairs(pairs), 2);
    EXPECT_TRUE(isPairInList(pairs, 0, 1));
    EXPECT_TRUE(isPairInList(pairs, 1, 0));
}

TEST_F(TestNeighborList, OneDirection) {
    // maxcutoff is 1.2, system is 4.8 x 5 x 5.1
    auto &ctx = kernel->getKernelContext();
    ctx.setBoxSize(1.2, 1.1, 2.8);
    ctx.setPeriodicBoundary(false, false, true);
    readdy::kernel::cpu::util::Config conf;
    cpum::NeighborList list(&ctx, &conf);
    list.setupBoxes();
    // Add three particles, one of which is in the neighborhood of the other two
    const auto particles = std::vector<m::Particle>{
            m::Particle(0, 0, -1.1, typeIdA), m::Particle(0, 0, .4, typeIdA), m::Particle(0, 0, 1.1, typeIdA)
    };
    scpum::ParticleData data;
    data.addParticles(particles);
    list.fillBoxes(data);
    auto pairs = list.pairs.get();
    EXPECT_EQ(getNumberPairs(pairs), 4);
    EXPECT_TRUE(isPairInList(pairs, 0, 2));
    EXPECT_TRUE(isPairInList(pairs, 2, 0));
    EXPECT_TRUE(isPairInList(pairs, 1, 2));
    EXPECT_TRUE(isPairInList(pairs, 2, 1));
    EXPECT_FALSE(isPairInList(pairs, 0, 1));
    EXPECT_FALSE(isPairInList(pairs, 1, 0));
}

TEST_F(TestNeighborList, AllNeighborsInCutoffSphere) {
    // maxcutoff is 1.2, system is 4 x 4 x 4, all directions periodic
    auto &ctx = kernel->getKernelContext();
    ctx.setBoxSize(4, 4, 4);
    ctx.setPeriodicBoundary(true, true, true);
    readdy::kernel::cpu::util::Config conf;
    cpum::NeighborList list(&ctx, &conf);
    list.setupBoxes();
    // Create a few particles. In this box setup, all particles are neighbors.
    const auto particles = std::vector<m::Particle>{
            m::Particle(0, 0, 0, typeIdA), m::Particle(0, 0, 0, typeIdA), m::Particle(.3, 0, 0, typeIdA),
            m::Particle(0, .3, -.3, typeIdA), m::Particle(-.3, 0, .3, typeIdA), m::Particle(.3, -.3, 0, typeIdA)
    };
    scpum::ParticleData data;
    data.addParticles(particles);
    list.fillBoxes(data);
    auto pairs = list.pairs.get();
    EXPECT_EQ(getNumberPairs(pairs), 30);
    for (size_t i = 0; i < 6; ++i) {
        for (size_t j = i + 1; j < 6; ++j) {
            EXPECT_TRUE(isPairInList(pairs, i, j)) << "Particles " << i << " and " << j << " were not neighbors.";
            EXPECT_TRUE(isPairInList(pairs, j, i)) << "Particles " << j << " and " << i << " were not neighbors.";
        }
    }
}
}