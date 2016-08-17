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
#include <readdy/testing/KernelContextMock.h>

using namespace readdy::kernel::singlecpu::model;

TEST(NeighborList, Naive) {
    unsigned int n_particles = 20;
    SingleCPUParticleData data;
    for (unsigned int i = 0; i < n_particles; ++i) {
        data.addParticle({(double) i, (double) i, (double) i, 5});
    }
    NaiveSingleCPUNeighborList list;
    list.create(data);
    EXPECT_EQ(((n_particles - 1) * n_particles) / 2, std::distance(list.begin(), list.end()));

    for(auto i = 0; i < n_particles; ++i) {
        for(auto j = i+1; j < n_particles; ++j) {
            EXPECT_TRUE(std::find(list.begin(), list.end(), ParticleIndexPair(j,i)) != list.end());
        }
    }
 }

TEST(NeighborList, NotThatNaive) {
    // Check very small system that only fits one box with many particles
    // Periodic and not periodic
    // Test that potentials are considered correctly w.r.t. cutoff via context mock
    readdy::testing::KernelContextMock ctx;
    readdy::testing::KernelContextMock const * const ctxPtr = &ctx;
    NotThatNaiveSingleCPUNeighborList<std::unordered_set<ParticleIndexPair, ParticleIndexPairHasher>,readdy::testing::KernelContextMock> list(ctxPtr);
}
