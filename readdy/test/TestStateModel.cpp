/**
 * @file TestStateModel.cpp
 * @brief Test the methods that manipulate time-dependent simulation data for kernels.
 * @author chrisfroe
 * @date 25.08.16
 * @todo check force calculation through periodic boundary
 */
#include <gtest/gtest.h>
#include <readdy/testing/Utils.h>
#include <readdy/testing/KernelTest.h>

namespace m = readdy::model;

namespace {

class TestStateModel : public KernelTest {

};

TEST_P(TestStateModel, CalculateForcesTwoParticles) {
    m::KernelContext &ctx = kernel->getKernelContext();
    auto &stateModel = kernel->getKernelStateModel();

    auto obs = kernel->createAndConnectObservable<m::ForcesObservable>(1);
    // two A particles with radius 1. -> cutoff 2, distance 1.8 -> r-r_0 = 0.2 -> force = 0.2
    ctx.setDiffusionConstant("A", 1.0);
    ctx.setParticleRadius("A", 1.0);
    ctx.setBoxSize(4., 4., 4.);
    ctx.setPeriodicBoundary(false, false, false);

    auto pot = kernel->createPotentialAs<m::potentials::HarmonicRepulsion>();
    pot->setForceConstant(1.);
    ctx.registerPotential(pot.get(), "A", "A");

    ctx.configure();
    auto typeIdA = ctx.getParticleTypeID("A");
    auto twoParticles = std::vector<m::Particle> {m::Particle(0., 0., 0., typeIdA), m::Particle(0., 0., 1.8, typeIdA)};
    stateModel.addParticles(twoParticles);
    stateModel.updateNeighborList();
    stateModel.calculateForces();
    stateModel.calculateForces(); // calculating twice should yield the same result. force and energy must not accumulate
    // check results
    std::get<0>(obs)->evaluate();
    auto forcesIt = std::get<0>(obs)->getResult().begin();
    EXPECT_VEC3_EQ(*forcesIt, m::Vec3(0, 0, -0.2));
    ++forcesIt;
    EXPECT_VEC3_EQ(*forcesIt, m::Vec3(0, 0, 0.2));
    EXPECT_DOUBLE_EQ(stateModel.getEnergy(), 0.02);
}

TEST_P(TestStateModel, CalculateForcesRepulsion) {
    m::KernelContext &ctx = kernel->getKernelContext();
    auto &stateModel = kernel->getKernelStateModel();

    // similar situation as before but now with repulsion between A and B
    auto obs = kernel->createAndConnectObservable<m::ForcesObservable>(1);
    ctx.setDiffusionConstant("A", 1.0);
    ctx.setDiffusionConstant("B", 1.0);
    ctx.setParticleRadius("A", 1.0);
    ctx.setParticleRadius("B", 2.0);
    ctx.setBoxSize(10., 10., 10.);
    ctx.setPeriodicBoundary(true, true, false);
    auto pot = kernel->createPotentialAs<m::potentials::HarmonicRepulsion>();
    {
        pot->setForceConstant(1.);
        ctx.registerPotential(pot.get(), "A", "B");
    }
    ctx.configure();
    auto typeIdA = ctx.getParticleTypeID("A");
    auto typeIdB = ctx.getParticleTypeID("B");
    /**
     * There are 6 particles. 0-2 are A particles. 3-5 are B particles.
     * The second B particle is a bit further away
     * -> There are forces between all AB pairs except with the second B particle, namely particle 4
     */
    auto particlesA = std::vector<m::Particle> {
            m::Particle(0, 0, 0, typeIdA), m::Particle(0, 0.8, 0, typeIdA), m::Particle(0.2, 0, -0.2, typeIdA)
    };
    auto particlesB = std::vector<m::Particle> {
            m::Particle(0, 0, 0.1, typeIdB), m::Particle(-1.8, -4.0, 1.8, typeIdB), m::Particle(0.5, 0, 0, typeIdB)
    };
    stateModel.addParticles(particlesA);
    stateModel.addParticles(particlesB);
    stateModel.updateNeighborList();
    stateModel.calculateForces();
    // handcalculated expectations
    const double energy03 = 4.205;
    const double energy05 = 3.125;
    const double energy13 = 2.4063226755104354;
    const double energy15 = 2.1148056603830194;
    const double energy23 = 3.4833346173608031;
    const double energy25 = 3.4833346173608031;
    const m::Vec3 force03(0, 0, -2.9);
    const m::Vec3 force05(-2.5, 0, 0);
    const m::Vec3 force13(0, 2.1768336301410027, -0.27210420376762534);
    const m::Vec3 force15(-1.08999682000954, 1.743994912015264, 0);
    const m::Vec3 force23(1.4641005886756873, 0, -2.1961508830135306);
    const m::Vec3 force25(-2.1961508830135306, 0, -1.4641005886756873);
    // check results
    std::get<0>(obs)->evaluate();
    auto forcesIt = std::get<0>(obs)->getResult().begin();
    EXPECT_VEC3_EQ(*forcesIt, force03 + force05) << "force on particle 0 = force03 + force05";
    ++forcesIt;
    EXPECT_VEC3_EQ(*forcesIt, force13 + force15) << "force on particle 1 = force13 + force15";
    ++forcesIt;
    EXPECT_VEC3_EQ(*forcesIt, force23 + force25) << "force on particle 2 = force23 + force25";
    ++forcesIt;
    EXPECT_VEC3_EQ(*forcesIt, (-1. * force03) - force13 - force23)
                        << "force on particle 3 = - force03 - force13 - force23";
    ++forcesIt;
    EXPECT_VEC3_EQ(*forcesIt, m::Vec3(0, 0, 0)) << "force on particle 4 = 0";
    ++forcesIt;
    EXPECT_VEC3_EQ(*forcesIt, (-1. * force05) - force15 - force25)
                        << "force on particle 5 = - force05 - force15 - force25";
    const double totalEnergy = energy03 + energy05 + energy13 + energy15 + energy23 + energy25;
    EXPECT_DOUBLE_EQ(stateModel.getEnergy(), totalEnergy);
}

TEST_P(TestStateModel, CalculateForcesNoForces) {
    m::KernelContext &ctx = kernel->getKernelContext();
    auto &stateModel = kernel->getKernelStateModel();

    // several particles without potentials -> forces must all be zero
    auto obs = kernel->createAndConnectObservable<m::ForcesObservable>(1);
    ctx.setDiffusionConstant("A", 1.0);
    ctx.setDiffusionConstant("B", 1.0);
    ctx.setParticleRadius("A", 1.0);
    ctx.setParticleRadius("B", 2.0);
    ctx.setBoxSize(4., 4., 4.);
    ctx.setPeriodicBoundary(false, false, false);
    ctx.configure();
    auto typeIdA = ctx.getParticleTypeID("A");
    auto typeIdB = ctx.getParticleTypeID("B");
    auto particlesA = std::vector<m::Particle> {
            m::Particle(0, 0, 0, typeIdA), m::Particle(0, 0.8, 0, typeIdA), m::Particle(0.2, 0, -0.2, typeIdA)
    };
    auto particlesB = std::vector<m::Particle> {
            m::Particle(0, 0, 0.1, typeIdB), m::Particle(-1.8, 0, 1.8, typeIdB), m::Particle(0.5, 0, 0, typeIdB)
    };
    stateModel.addParticles(particlesA);
    stateModel.addParticles(particlesB);
    stateModel.updateNeighborList();
    stateModel.calculateForces();
    // check results
    std::get<0>(obs)->evaluate();
    for (auto &&force : std::get<0>(obs)->getResult()) {
        EXPECT_VEC3_EQ(force, m::Vec3(0, 0, 0));
    }
}

INSTANTIATE_TEST_CASE_P(TestStateModel, TestStateModel,
                        ::testing::ValuesIn(readdy::testing::getKernelsToTest()));
}