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
            auto& ctx = kernel->getKernelContext();
            ctx.setDiffusionConstant("A", 1.0);
            double eductDistance = 1.2;
            ctx.registerReaction(kernel->createFusionReaction("test", "A", "A", "A", 0., eductDistance));
            readdy::testing::NOOPPotentialOrder2 pot(1.1, 0., 0.);
            ctx.registerOrder2Potential(&pot, "A", "A");
            typeIdA = ctx.getParticleTypeID("A");
            ctx.configure();
        }

    };

    auto isPairInList = [](std::unordered_map<unsigned long, std::vector<unsigned long>> *pairs, unsigned long idx1, unsigned long idx2) {
        return (std::find(pairs->operator[](idx1).cbegin(), pairs->operator[](idx1).cend(), idx2) != pairs->operator[](idx1).cend());
    };

    auto getNumberPairs = [](std::unordered_map<unsigned long, std::vector<unsigned long>> *pairs) {
        auto pairsIt = pairs->cbegin();
        size_t numberPairs = 0;
        while (pairsIt != pairs->cend()) {
            numberPairs += pairsIt->second.size();
            ++pairsIt;
        }
        return numberPairs;
    };

    TEST_F(TestNeighborList, ThreeBoxesPeriodicAxis) {
        // maxcutoff is 1.2 , system is 3.6 x 2 x 2, i.e. there are three boxes along the periodic axis
        auto& ctx = kernel->getKernelContext();
        ctx.setBoxSize(3.6, 2, 2);
        ctx.setPeriodicBoundary(true, false, false);
        cpum::NeighborList list(&ctx);
        list.setupBoxes();
        auto boxes = list.getBoxes();
        EXPECT_EQ(boxes.size(), 3);
        for (size_t b = 0; b < 3; ++b) {
            EXPECT_EQ(boxes[b].j, 0);
            EXPECT_EQ(boxes[b].k, 0);
            EXPECT_EQ(boxes[b].neighboringBoxes.size(), 2);
            EXPECT_TRUE((*boxes[b].neighboringBoxes[0]).id != boxes[b].id) << "A Box should not have itself as a neighbor.";
        }
        // now create three particles. The resulting neighborlist should contain six pairs
        const auto threeParticles = std::vector<m::Particle>{
                m::Particle(0, 0, 0, typeIdA), m::Particle(0, 0, 0, typeIdA), m::Particle(1.6, 0, 0, typeIdA)};
        scpum::SingleCPUParticleData data;
        data.addParticles(threeParticles);
        list.fillBoxes(data);
        // check result
        auto pairs = list.pairs.get();
        EXPECT_EQ(getNumberPairs(pairs), 6);
        EXPECT_TRUE(isPairInList(pairs, 0, 1));
        EXPECT_TRUE(isPairInList(pairs, 1, 0));
        EXPECT_TRUE(isPairInList(pairs, 0, 2));
        EXPECT_TRUE(isPairInList(pairs, 2, 0));
        EXPECT_TRUE(isPairInList(pairs, 1, 2));
        EXPECT_TRUE(isPairInList(pairs, 2, 1));
    }

    TEST_F(TestNeighborList, ThreeBoxesNonPeriodic) {
        // maxcutoff is 1.2, system is 1.5 x 4 x 1.5, non-periodic, three boxes
        auto& ctx = kernel->getKernelContext();
        ctx.setBoxSize(1.5, 4, 1.5);
        ctx.setPeriodicBoundary(false, false, false);
        cpum::NeighborList list(&ctx);
        list.setupBoxes();
        auto boxes = list.getBoxes();
        EXPECT_EQ(boxes.size(), 3);
        // Add three particles, two are in one outer box, the third on the other end and thus no neighbor
        const auto particles = std::vector<m::Particle>{
                m::Particle(0, -1.8, 0, typeIdA), m::Particle(0, -1.8, 0, typeIdA), m::Particle(0, 1.8, 0, typeIdA)
        };
        scpum::SingleCPUParticleData data;
        data.addParticles(particles);
        list.fillBoxes(data);
        auto pairs = list.pairs.get();
        EXPECT_EQ(getNumberPairs(pairs), 2);
        EXPECT_TRUE(isPairInList(pairs, 0, 1));
        EXPECT_TRUE(isPairInList(pairs, 1, 0));
    }

    TEST_F(TestNeighborList, 64BoxesAllPeriodic) {
        // maxcutoff is 1.2, system is 4.8 x 5 x 5.1, all periodic, i.e. 64 boxes each with 26 neighbors
        auto& ctx = kernel->getKernelContext();
        ctx.setBoxSize(4.8, 5, 5.1);
        ctx.setPeriodicBoundary(true, true, true);
        cpum::NeighborList list(&ctx);
        list.setupBoxes();
        auto boxes = list.getBoxes();
        EXPECT_EQ(boxes.size(), 64);
        for (auto &&box : boxes) {
            EXPECT_EQ(box.neighboringBoxes.size(), 26);
        }
        // Add three particles, one of which is in the neighborhood of the other two
        const auto particles = std::vector<m::Particle>{
                m::Particle(-2.1, -2.4, -2.4, typeIdA), m::Particle(1, 1, 1, typeIdA), m::Particle(2.1, 2.4, 2.4, typeIdA)
        };
        scpum::SingleCPUParticleData data;
        data.addParticles(particles);
        list.fillBoxes(data);
        auto pairs = list.pairs.get();
        EXPECT_EQ(getNumberPairs(pairs), 4);
        EXPECT_TRUE(isPairInList(pairs, 0, 2));
        EXPECT_TRUE(isPairInList(pairs, 2, 0));
        EXPECT_TRUE(isPairInList(pairs, 1, 2));
        EXPECT_TRUE(isPairInList(pairs, 2, 1));
    }

    TEST_F(TestNeighborList, 27BoxesAllPeriodic) {
        // maxcutoff is 1.2, system is 4 x 4 x 4, all directions periodic, i.e. 27 boxes each with 26 neighbors
        auto& ctx = kernel->getKernelContext();
        ctx.setBoxSize(4, 4, 4);
        ctx.setPeriodicBoundary(true, true, true);
        cpum::NeighborList list(&ctx);
        list.setupBoxes();
        auto boxes = list.getBoxes();
        EXPECT_EQ(boxes.size(), 27);
        for (auto &&box :boxes) {
            EXPECT_EQ(box.neighboringBoxes.size(), 26);
        }
        // Create a few particles. In this box setup, all particles are neighbors.
        const auto particles = std::vector<m::Particle>{
                m::Particle(0, 0, 0, typeIdA), m::Particle(0, 0, 0, typeIdA), m::Particle(1.6, 0, 0, typeIdA),
                m::Particle(0, 1.6, -1.6, typeIdA), m::Particle(-1.6, 0, 1.6, typeIdA), m::Particle(1.6, -1.6, 0, typeIdA)
        };
        scpum::SingleCPUParticleData data;
        data.addParticles(particles);
        list.fillBoxes(data);
        auto pairs = list.pairs.get();
        EXPECT_EQ(getNumberPairs(pairs), 30);
        for (size_t i = 0; i < 6; ++i) {
            for (size_t j = i + 1; j < 6; ++j) {
                EXPECT_TRUE(isPairInList(pairs, i, j));
                EXPECT_TRUE(isPairInList(pairs, j, i));
            }
        }
    }
}