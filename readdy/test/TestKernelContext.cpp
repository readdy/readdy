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
#include <readdy/model/potentials/PotentialsOrder1.h>
#include <readdy/model/potentials/PotentialsOrder2.h>
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
    auto p1 = std::make_unique<readdy::testing::NOOPPotentialOrder2>();
    ctx.registerOrder2Potential(p1.get(), "a", "b");
    ctx.registerOrder2Potential(p1.get(), "b", "a");
    ctx.configure();
    auto vector = ctx.getOrder2Potentials("b", "a");
    EXPECT_EQ(vector.size(), 2);
}

TEST_P(TestKernelContextWithKernels, PotentialOrder1Map) {
    auto kernel = readdy::plugin::KernelProvider::getInstance().create("SingleCPU");

    namespace rmp = readdy::model::potentials;

    kernel->getKernelContext().setDiffusionConstant("A", 1.0);
    kernel->getKernelContext().setDiffusionConstant("B", 3.0);
    kernel->getKernelContext().setDiffusionConstant("C", 4.0);
    kernel->getKernelContext().setDiffusionConstant("D", 2.0);

    kernel->getKernelContext().setParticleRadius("A", 1.0);
    kernel->getKernelContext().setParticleRadius("B", 2.0);
    kernel->getKernelContext().setParticleRadius("C", 3.0);
    kernel->getKernelContext().setParticleRadius("D", 4.0);

    boost::uuids::uuid uuid1_1, uuid1_11;
    boost::uuids::uuid uuid2_1, uuid2_2;
    {
        auto pot1 = kernel->createPotentialAs<rmp::CubePotential>();
        auto pot11 = kernel->createPotentialAs<rmp::CubePotential>();
        auto pot2 = kernel->createPotentialAs<rmp::HarmonicRepulsion>();
        uuid1_1 = pot1->getId();

        kernel->getKernelContext().registerOrder1Potential(pot1.get(), "A");
        kernel->getKernelContext().registerOrder1Potential(pot1.get(), "C");
        kernel->getKernelContext().registerOrder1Potential(pot1.get(), "D");
        kernel->getKernelContext().registerOrder1Potential(pot1.get(), "C");

        uuid1_11 = kernel->getKernelContext().registerOrder1Potential(pot11.get(), "B");

        uuid2_1 = kernel->getKernelContext().registerOrder2Potential(pot2.get(), "A", "C");
        uuid2_2 = kernel->getKernelContext().registerOrder2Potential(pot2.get(), "B", "C");
        kernel->getKernelContext().configure();
        EXPECT_EQ(pot2->getId(), uuid2_1);
        EXPECT_EQ(pot2->getId(), uuid2_2);
    }
    // test that order 1 potentials are set up correctly
    {
        {
            const auto &pot1_A = kernel->getKernelContext().getOrder1Potentials("A");
            EXPECT_EQ(pot1_A.size(), 1);
            EXPECT_EQ(pot1_A[0]->getId(), uuid1_1);
            EXPECT_EQ(pot1_A[0]->getName(), rmp::getPotentialName<rmp::CubePotential>());
            EXPECT_EQ(dynamic_cast<rmp::CubePotential *>(pot1_A[0])->getParticleRadius(), 1.0);
        }
        {
            const auto &pot1_B = kernel->getKernelContext().getOrder1Potentials("B");
            EXPECT_EQ(pot1_B.size(), 1);
            EXPECT_EQ(pot1_B[0]->getId(), uuid1_11);
            EXPECT_EQ(pot1_B[0]->getName(), rmp::getPotentialName<rmp::CubePotential>());
            EXPECT_EQ(dynamic_cast<rmp::CubePotential *>(pot1_B[0])->getParticleRadius(), 2.0);
        }
        {
            const auto &pot1_C = kernel->getKernelContext().getOrder1Potentials("C");
            EXPECT_EQ(pot1_C.size(), 2);
            for (auto &&ptr : pot1_C) {
                EXPECT_EQ(ptr->getId(), uuid1_1);
                EXPECT_EQ(ptr->getName(), rmp::getPotentialName<rmp::CubePotential>());
                EXPECT_EQ(dynamic_cast<rmp::CubePotential *>(ptr)->getParticleRadius(), 3.0);
            }
        }
        {
            const auto &pot1_D = kernel->getKernelContext().getOrder1Potentials("D");
            EXPECT_EQ(pot1_D.size(), 1);
            EXPECT_EQ(pot1_D[0]->getId(), uuid1_1);
            EXPECT_EQ(pot1_D[0]->getName(), rmp::getPotentialName<rmp::CubePotential>());
            EXPECT_EQ(dynamic_cast<rmp::CubePotential *>(pot1_D[0])->getParticleRadius(), 4.0);
        }
    }
    // test that order 2 potentials are set up correctly
    {
        EXPECT_EQ(kernel->getKernelContext().getOrder2Potentials("A", "A").size(), 0);
        EXPECT_EQ(kernel->getKernelContext().getOrder2Potentials("A", "B").size(), 0);
        EXPECT_EQ(kernel->getKernelContext().getOrder2Potentials("A", "D").size(), 0);
        EXPECT_EQ(kernel->getKernelContext().getOrder2Potentials("B", "B").size(), 0);
        EXPECT_EQ(kernel->getKernelContext().getOrder2Potentials("B", "D").size(), 0);

        EXPECT_EQ(kernel->getKernelContext().getOrder2Potentials("A", "C").size(), 1);
        EXPECT_EQ(kernel->getKernelContext().getOrder2Potentials("B", "C").size(), 1);
        {
            const auto &pot2_AC = kernel->getKernelContext().getOrder2Potentials("A", "C");
            EXPECT_EQ(pot2_AC[0]->getId(), uuid2_1);
            EXPECT_EQ(pot2_AC[0]->getName(), rmp::getPotentialName<rmp::HarmonicRepulsion>());
            EXPECT_EQ(dynamic_cast<rmp::HarmonicRepulsion *>(pot2_AC[0])->getSumOfParticleRadii(), 1 + 3);
        }
        {
            const auto &pot2_BC = kernel->getKernelContext().getOrder2Potentials("B", "C");
            EXPECT_EQ(pot2_BC[0]->getId(), uuid2_2);
            EXPECT_EQ(pot2_BC[0]->getName(), rmp::getPotentialName<rmp::HarmonicRepulsion>());
            EXPECT_EQ(dynamic_cast<rmp::HarmonicRepulsion *>(pot2_BC[0])->getSumOfParticleRadii(), 2 + 3);
        }
    }

    // now remove
    {
        // only one potential for particle type B has a different uuid
        kernel->getKernelContext().deregisterPotential(uuid1_1);
        kernel->getKernelContext().configure();
        EXPECT_EQ(kernel->getKernelContext().getOrder1Potentials("A").size(), 0);
        EXPECT_EQ(kernel->getKernelContext().getOrder1Potentials("B").size(), 1);
        EXPECT_EQ(kernel->getKernelContext().getOrder1Potentials("C").size(), 0);
        EXPECT_EQ(kernel->getKernelContext().getOrder1Potentials("D").size(), 0);

        // both potentials have same uuid, so they should be discarded both
        kernel->getKernelContext().deregisterPotential(uuid2_2);
        kernel->getKernelContext().configure();
        EXPECT_EQ(kernel->getKernelContext().getOrder2Potentials("A", "C").size(), 0);
        EXPECT_EQ(kernel->getKernelContext().getOrder2Potentials("B", "C").size(), 0);
    }

}

INSTANTIATE_TEST_CASE_P(TestKernelContext, TestKernelContextWithKernels,
                        ::testing::ValuesIn(readdy::testing::getKernelsToTest()));

}