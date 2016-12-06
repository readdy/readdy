/**
 * Test different contexts w.r.t. boxsize and periodicity, perform setupBoxes() and see if that worked.
 * Then add a few particles and perform fillBoxes(). The result of this is the actual neighborlist, which is
 * checked as well.
 *
 * @file SingleCPUTestNeighborList.cpp
 * @brief Single cpu neighbor list tests
 * @author clonker
 * @author chrisfroe
 * @date 09.06.16
 */

#include <gtest/gtest.h>
#include <readdy/kernel/singlecpu/model/SCPUNeighborList.h>
#include <readdy/kernel/singlecpu/SCPUKernel.h>
#include <readdy/testing/NOOPPotential.h>

namespace scpu = readdy::kernel::scpu;
namespace scpum = scpu::model;
namespace m = readdy::model;

using kernel_t = readdy::kernel::scpu::SCPUKernel;

struct NeighborListTest : ::testing::Test {

    std::unique_ptr<kernel_t> kernel;
    unsigned int typeIdA;

    NeighborListTest() : kernel(std::make_unique<kernel_t>()) {
        readdy::model::KernelContext &ctx = kernel->getKernelContext();
        ctx.setDiffusionConstant("A", 1.0);
        double eductDistance = 1.2;
        kernel->registerReaction<readdy::model::reactions::Fusion>("test", "A", "A", "A", 0., eductDistance);
        ctx.registerPotential(std::make_unique<readdy::testing::NOOPPotentialOrder2>(1.1,0,0), "A", "A");
        typeIdA = ctx.getParticleTypeID("A");
        ctx.configure();
    }
};

TEST(NeighborList, Naive) {
    unsigned int n_particles = 20;
    scpum::SCPUParticleData data;
    for (unsigned int i = 0; i < n_particles; ++i) {
        data.addParticle({(double) i, (double) i, (double) i, 5});
    }
    scpum::SCPUNaiveNeighborList list;
    list.create(data);
    EXPECT_EQ(((n_particles - 1) * n_particles) / 2, std::distance(list.begin(), list.end()));

    for (auto i = 0; i < n_particles; ++i) {
        for (auto j = i + 1; j < n_particles; ++j) {
            EXPECT_TRUE(std::find(list.begin(), list.end(), scpum::ParticleIndexPair(j, i)) != list.end());
        }
    }
}

TEST_F(NeighborListTest, ThreeBoxesPeriodicAxis) {
    // maxcutoff is 1.2 , system is 3.6 x 2 x 2, i.e. there are three cells along the periodic axis
    auto &ctx = kernel->getKernelContext();
    ctx.setBoxSize(3.6, 2, 2);
    ctx.setPeriodicBoundary(true, false, false);
    scpum::SCPUNotThatNaiveNeighborList<> list(&ctx);
    list.setupBoxes();
    auto boxes = list.getBoxes();
    EXPECT_EQ(boxes.size(), 3);
    for (size_t b = 0; b < 3; ++b) {
        EXPECT_EQ(boxes[b].j, 0);
        EXPECT_EQ(boxes[b].k, 0);
        EXPECT_EQ(boxes[b].neighbors.size(), 1);
        EXPECT_TRUE(*boxes[b].neighbors[0] != boxes[b])
                            << "This compares ids. A Box should not have itself as a neighbor.";
    }
    // now create three particles. The resulting neighborlist should contain three pairs
    const auto threeParticles = std::vector<m::Particle>{
            m::Particle(0, 0, 0, typeIdA), m::Particle(0, 0, 0, typeIdA), m::Particle(1.6, 0, 0, typeIdA)};
    scpum::SCPUParticleData data;
    data.addParticles(threeParticles);
    list.fillBoxes(data);
    EXPECT_EQ(std::distance(list.cbegin(), list.cend()), 3);
    EXPECT_TRUE(std::find(list.cbegin(), list.cend(), scpum::ParticleIndexPair(0, 1)) != list.cend())
                        << "neighborlist should contain (0,1)";
    EXPECT_TRUE(std::find(list.cbegin(), list.cend(), scpum::ParticleIndexPair(0, 2)) != list.cend())
                        << "neighborlist should contain (0,2)";
    EXPECT_TRUE(std::find(list.cbegin(), list.cend(), scpum::ParticleIndexPair(1, 2)) != list.cend())
                        << "neighborlist should contain (1,2)";


}

TEST_F(NeighborListTest, 27BoxesAllPeriodic) {
    // maxcutoff is 1.2, system is 4 x 4 x 4, all directions periodic, i.e. 27 cells each with 13 neighbors
    auto &ctx = kernel->getKernelContext();
    ctx.setBoxSize(4, 4, 4);
    ctx.setPeriodicBoundary(true, true, true);
    scpum::SCPUNotThatNaiveNeighborList<> list(&ctx);
    list.setupBoxes();
    auto boxes = list.getBoxes();
    EXPECT_EQ(boxes.size(), 27);
    for (auto &&box : boxes) {
        EXPECT_EQ(box.neighbors.size(), 13);
    }
    // Create a few particles. In this box setup, all particles are neighbors.
    const auto particles = std::vector<m::Particle>{
            m::Particle(0, 0, 0, typeIdA), m::Particle(0, 0, 0, typeIdA), m::Particle(1.6, 0, 0, typeIdA),
            m::Particle(0, 1.6, -1.6, typeIdA), m::Particle(-1.6, 0, 1.6, typeIdA), m::Particle(1.6, -1.6, 0, typeIdA)
    };
    scpum::SCPUParticleData data;
    data.addParticles(particles);
    list.fillBoxes(data);
    EXPECT_EQ(std::distance(list.cbegin(), list.cend()), 15);
    for (size_t i = 0; i < 6; ++i) {
        for (size_t j = i + 1; j < 6; ++j) {
            EXPECT_TRUE(std::find(list.cbegin(), list.cend(), scpum::ParticleIndexPair(i, j)) != list.cend());
        }
    }
}

TEST_F(NeighborListTest, 64BoxesAllPeriodic) {
    // maxcutoff is 1.2, system is 4.8 x 5 x 5.1, all periodic, i.e. 64 cells each with 13 neighbors
    auto &ctx = kernel->getKernelContext();
    ctx.setBoxSize(4.8, 5, 5.1);
    ctx.setPeriodicBoundary(true, true, true);
    scpum::SCPUNotThatNaiveNeighborList<> list(&ctx);
    list.setupBoxes();
    auto boxes = list.getBoxes();
    EXPECT_EQ(boxes.size(), 64);
    for (auto &&box : boxes) {
        EXPECT_EQ(box.neighbors.size(), 13);
    }
    // Add three particles, one of which is in the neighborhood of the other two
    const auto particles = std::vector<m::Particle>{
            m::Particle(-2.1, -2.4, -2.4, typeIdA), m::Particle(1, 1, 1, typeIdA), m::Particle(2.1, 2.4, 2.4, typeIdA)
    };
    scpum::SCPUParticleData data;
    data.addParticles(particles);
    list.fillBoxes(data);
    EXPECT_EQ(std::distance(list.cbegin(), list.cend()), 2);
    EXPECT_TRUE(std::find(list.cbegin(), list.cend(), scpum::ParticleIndexPair(0, 2)) != list.cend())
                        << "0 is neighbor of 2";
    EXPECT_TRUE(std::find(list.cbegin(), list.cend(), scpum::ParticleIndexPair(1, 2)) != list.cend())
                        << "1 is neighbor of 2";
}


TEST_F(NeighborListTest, ThreeBoxesNonPeriodic) {
    // maxcutoff is 1.2, system is 1.5 x 4 x 1.5, non-periodic, three cells
    auto &ctx = kernel->getKernelContext();
    ctx.setBoxSize(1.5, 4, 1.5);
    ctx.setPeriodicBoundary(false, false, false);
    scpum::SCPUNotThatNaiveNeighborList<> list(&ctx);
    list.setupBoxes();
    auto boxes = list.getBoxes();
    EXPECT_EQ(boxes.size(), 3);
    // Add three particles, two are in one outer box, the third on the other end and thus no neighbor
    const auto particles = std::vector<m::Particle>{
            m::Particle(0, -1.8, 0, typeIdA), m::Particle(0, -1.8, 0, typeIdA), m::Particle(0, 1.8, 0, typeIdA)
    };
    scpum::SCPUParticleData data;
    data.addParticles(particles);
    list.fillBoxes(data);
    EXPECT_EQ(std::distance(list.cbegin(), list.cend()), 1);
    EXPECT_TRUE(std::find(list.cbegin(), list.cend(), scpum::ParticleIndexPair(0, 1)) != list.cend())
                        << "only 0 and 1 are neighbors";
}


