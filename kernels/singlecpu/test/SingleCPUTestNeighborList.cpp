/**
 * << detailed description >>
 *
 * @file SingleCPUTestNeighborList.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 09.06.16
 */

#include <gtest/gtest.h>
#include <readdy/kernel/singlecpu/model/SingleCPUNeighborList.h>
#include <readdy/kernel/singlecpu/SingleCPUKernel.h>
#include <readdy/kernel/singlecpu/reactions/SingleCPUReactionFactory.h>
#include <readdy/testing/NOOPPotential.h>

namespace scpu = readdy::kernel::singlecpu;
namespace scpum = scpu::model;

TEST(NeighborList, Naive) {
    unsigned int n_particles = 20;
    scpum::SingleCPUParticleData data;
    for (unsigned int i = 0; i < n_particles; ++i) {
        data.addParticle({(double) i, (double) i, (double) i, 5});
    }
    scpum::NaiveSingleCPUNeighborList list;
    list.create(data);
    EXPECT_EQ(((n_particles - 1) * n_particles) / 2, std::distance(list.begin(), list.end()));

    for(auto i = 0; i < n_particles; ++i) {
        for(auto j = i+1; j < n_particles; ++j) {
            EXPECT_TRUE(std::find(list.begin(), list.end(), scpum::ParticleIndexPair(j,i)) != list.end());
        }
    }
 }

TEST(NeighborList, NotThatNaive) {
    // Check very small system that only fits one box with many particles
    // Periodic and not periodic
    // Test that potentials are considered correctly w.r.t. cutoff via context mock
    // First check setupBoxes
    scpu::SingleCPUKernel kernel;
    scpu::reactions::SingleCPUReactionFactory reactionFactory(&kernel);
    readdy::model::KernelContext ctx(&reactionFactory);
    ctx.setDiffusionConstant("A", 1.0);
    double eductDistance = 5.0;
    ctx.registerFusionReaction("test", "A", "A", "A", 1.0, eductDistance);
    const readdy::testing::NOOPPotential pot;
    ctx.registerOrder2Potential(&pot, "A", "A");
    scpum::NotThatNaiveSingleCPUNeighborList<> list(&ctx);
}
