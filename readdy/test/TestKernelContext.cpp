/**
 * Testfile containing tests for the KernelContext class.
 *
 * @file TestKernelContext.cpp
 * @brief Testfile for KernelContext.
 * @author clonker
 * @date 19.04.16
 */


#include <readdy/common/Utils.h>
#include <readdy/common/make_unique.h>
#include <readdy/model/Kernel.h>
#include <boost/algorithm/string.hpp>
#include <readdy/plugin/KernelProvider.h>
#include <readdy/testing/KernelTest.h>
#include <readdy/testing/Utils.h>
#include <readdy/testing/NOOPPotential.h>

namespace m = readdy::model;

namespace {

class TestKernelContext : public ::testing::Test {
protected:
    TestKernelContext() {}

    std::unique_ptr<readdy::model::reactions::ReactionFactory> reactionFactory{
            new readdy::model::reactions::ReactionFactory()};
};

class TestKernelContextWithKernels : public KernelTest {

};

TEST_F(TestKernelContext, SetGetKBT) {
    m::KernelContext ctx;
    ctx.setKBT(42);
    EXPECT_EQ(42, ctx.getKBT());
}

TEST_F(TestKernelContext, PeriodicBoundary) {
    m::KernelContext ctx;
    ctx.setPeriodicBoundary(true, false, true);
    auto boundary = ctx.getPeriodicBoundary();
    EXPECT_TRUE(boundary[0]);
    EXPECT_FALSE(boundary[1]);
    EXPECT_TRUE(boundary[2]);
}

TEST_F(TestKernelContext, BoxSize) {
    m::KernelContext ctx;
    ctx.setBoxSize(10, 11, 12);
    auto box_size = ctx.getBoxSize();
    EXPECT_EQ(box_size[0], 10);
    EXPECT_EQ(box_size[1], 11);
    EXPECT_EQ(box_size[2], 12);
}

TEST_F(TestKernelContext, PotentialOrder2Map) {
    m::KernelContext ctx;
    ctx.registerPotential(std::make_unique<readdy::testing::NOOPPotentialOrder2>(), "a", "b");
    ctx.registerPotential(std::make_unique<readdy::testing::NOOPPotentialOrder2>(), "b", "a");
    ctx.configure();
    auto vector = ctx.getOrder2Potentials("b", "a");
    EXPECT_EQ(vector.size(), 2);
}

TEST_P(TestKernelContextWithKernels, PotentialOrder1Map) {
    auto kernel = readdy::plugin::KernelProvider::getInstance().create("SingleCPU");

    namespace rmp = readdy::model::potentials;
    
    auto& ctx = kernel->getKernelContext();

    ctx.setDiffusionConstant("A", 1.0);
    ctx.setDiffusionConstant("B", 3.0);
    ctx.setDiffusionConstant("C", 4.0);
    ctx.setDiffusionConstant("D", 2.0);

    ctx.setParticleRadius("A", 1.0);
    ctx.setParticleRadius("B", 2.0);
    ctx.setParticleRadius("C", 3.0);
    ctx.setParticleRadius("D", 4.0);
    std::vector<short> idsToRemove;
    short uuid2_1, uuid2_2;
    {
        auto id1 = ctx.registerPotential(kernel->createPotentialAs<rmp::CubePotential>(), "A");
        auto id2 = ctx.registerPotential(kernel->createPotentialAs<rmp::CubePotential>(), "C");
        auto id3 = ctx.registerPotential(kernel->createPotentialAs<rmp::CubePotential>(), "D");
        auto id4 = ctx.registerPotential(kernel->createPotentialAs<rmp::CubePotential>(), "C");
        
        idsToRemove.push_back(id1);
        idsToRemove.push_back(id2);
        idsToRemove.push_back(id3);
        idsToRemove.push_back(id4);

        ctx.registerPotential(kernel->createPotentialAs<rmp::CubePotential>(), "B");

        uuid2_1 = ctx.registerPotential(kernel->createPotentialAs<rmp::HarmonicRepulsion>(), "A", "C");
        uuid2_2 = ctx.registerPotential(kernel->createPotentialAs<rmp::HarmonicRepulsion>(), "B", "C");
        ctx.configure();
    }
    // test that order 1 potentials are set up correctly
    {
        {
            const auto &pot1_A = ctx.getOrder1Potentials("A");
            EXPECT_EQ(pot1_A.size(), 1);
            EXPECT_EQ(pot1_A[0]->getName(), rmp::getPotentialName<rmp::CubePotential>());
            EXPECT_EQ(dynamic_cast<rmp::CubePotential *>(pot1_A[0])->getParticleRadius(), 1.0);
        }
        {
            const auto &pot1_B = ctx.getOrder1Potentials("B");
            EXPECT_EQ(pot1_B.size(), 1);
            EXPECT_EQ(pot1_B[0]->getName(), rmp::getPotentialName<rmp::CubePotential>());
            EXPECT_EQ(dynamic_cast<rmp::CubePotential *>(pot1_B[0])->getParticleRadius(), 2.0);
        }
        {
            const auto &pot1_C = ctx.getOrder1Potentials("C");
            EXPECT_EQ(pot1_C.size(), 2);
            for (auto &&ptr : pot1_C) {
                EXPECT_EQ(ptr->getName(), rmp::getPotentialName<rmp::CubePotential>());
                EXPECT_EQ(dynamic_cast<rmp::CubePotential *>(ptr)->getParticleRadius(), 3.0);
            }
        }
        {
            const auto &pot1_D = ctx.getOrder1Potentials("D");
            EXPECT_EQ(pot1_D.size(), 1);
            EXPECT_EQ(pot1_D[0]->getName(), rmp::getPotentialName<rmp::CubePotential>());
            EXPECT_EQ(dynamic_cast<rmp::CubePotential *>(pot1_D[0])->getParticleRadius(), 4.0);
        }
    }
    // test that order 2 potentials are set up correctly
    {
        EXPECT_EQ(ctx.getOrder2Potentials("A", "A").size(), 0);
        EXPECT_EQ(ctx.getOrder2Potentials("A", "B").size(), 0);
        EXPECT_EQ(ctx.getOrder2Potentials("A", "D").size(), 0);
        EXPECT_EQ(ctx.getOrder2Potentials("B", "B").size(), 0);
        EXPECT_EQ(ctx.getOrder2Potentials("B", "D").size(), 0);

        EXPECT_EQ(ctx.getOrder2Potentials("A", "C").size(), 1);
        EXPECT_EQ(ctx.getOrder2Potentials("B", "C").size(), 1);
        {
            const auto &pot2_AC = ctx.getOrder2Potentials("A", "C");
            EXPECT_EQ(pot2_AC[0]->getId(), uuid2_1);
            EXPECT_EQ(pot2_AC[0]->getName(), rmp::getPotentialName<rmp::HarmonicRepulsion>());
            EXPECT_EQ(dynamic_cast<rmp::HarmonicRepulsion *>(pot2_AC[0])->getSumOfParticleRadii(), 1 + 3);
        }
        {
            const auto &pot2_BC = ctx.getOrder2Potentials("B", "C");
            EXPECT_EQ(pot2_BC[0]->getId(), uuid2_2);
            EXPECT_EQ(pot2_BC[0]->getName(), rmp::getPotentialName<rmp::HarmonicRepulsion>());
            EXPECT_EQ(dynamic_cast<rmp::HarmonicRepulsion *>(pot2_BC[0])->getSumOfParticleRadii(), 2 + 3);
        }
    }

    // now remove
    std::for_each(idsToRemove.begin(), idsToRemove.end(), [&](const short id) {ctx.deregisterPotential(id);});
    {
        // only one potential for particle type B has a different uuid
        ctx.configure();
        EXPECT_EQ(ctx.getOrder1Potentials("A").size(), 0);
        EXPECT_EQ(ctx.getOrder1Potentials("B").size(), 1);
        EXPECT_EQ(ctx.getOrder1Potentials("C").size(), 0);
        EXPECT_EQ(ctx.getOrder1Potentials("D").size(), 0);

        // remove 2nd potential
        ctx.deregisterPotential(uuid2_2);
        ctx.configure();
        EXPECT_EQ(ctx.getOrder2Potentials("A", "C").size(), 1);
        EXPECT_EQ(ctx.getOrder2Potentials("B", "C").size(), 0);
    }

}

INSTANTIATE_TEST_CASE_P(TestKernelContext, TestKernelContextWithKernels,
                        ::testing::ValuesIn(readdy::testing::getKernelsToTest()));

}