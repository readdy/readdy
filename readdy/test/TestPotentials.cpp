/**
 * Test potentials of first and second order via small particle systems. First order potentials
 * often assure that particles stay in a certain volume of the simulation box, which is checked for.
 *
 * @file TestPotentials.cpp
 * @brief Check correct behavior of particles subject to potentials.
 * @author clonker
 * @author chrisfroe
 * @date 27.06.16
 */

#include <gtest/gtest.h>
#include <readdy/plugin/KernelProvider.h>
#include <readdy/testing/KernelTest.h>
#include <readdy/testing/Utils.h>

namespace {

class TestPotentials : public KernelTest {

};

void setupParticles(readdy::model::Kernel &kernel) {
    kernel.getKernelContext().setDiffusionConstant("A", 1);
    kernel.getKernelContext().setDiffusionConstant("B", .1);
    kernel.getKernelContext().setParticleRadius("A", .1);
    kernel.getKernelContext().setParticleRadius("B", .01);
    kernel.getKernelContext().setPeriodicBoundary(false, false, false);
    const unsigned int nParticlesA = 10;
    const unsigned int nParticlesB = 10;

    // add particles close to origin but not all exactly on 0,0,0
    for (auto i = 0; i < nParticlesA; ++i) kernel.addParticle("A", readdy::model::rnd::normal3());
    for (auto i = 0; i < nParticlesB; ++i) kernel.addParticle("B", readdy::model::rnd::normal3());
}

void run(readdy::model::Kernel &kernel) {
    unsigned int nSteps = 2000;
    auto &&integrator = kernel.createProgram<readdy::model::programs::EulerBDIntegrator>();
    auto &&nl = kernel.createProgram<readdy::model::programs::UpdateNeighborList>();
    auto &&forces = kernel.createProgram<readdy::model::programs::CalculateForces>();
    nl->execute();
    for (readdy::model::observables::time_step_type &&t = 0; t < nSteps; ++t) {
        forces->execute();
        integrator->execute();
        nl->execute();
    }
}

TEST_P(TestPotentials, TestParticlesStayInBox) {
    kernel->getKernelContext().setBoxSize(5, 5, 5);
    kernel->getKernelContext().setTimeStep(.005);

    setupParticles(*kernel);

    auto harmonicRepulsion = kernel->createPotentialAs<readdy::model::potentials::HarmonicRepulsion>();
    harmonicRepulsion->setForceConstant(.01);
    kernel->getKernelContext().registerPotential(std::move(harmonicRepulsion), "A", "B");
    std::array<std::string, 2> types{{"A", "B"}};
    for (auto t : types) {
        auto cubePot = kernel->createPotentialAs<readdy::model::potentials::CubePotential>();
        cubePot->setOrigin({-1, -1, -1});
        cubePot->setConsiderParticleRadius(true);
        cubePot->setExtent({2, 2, 2});
        cubePot->setForceConstant(10);
        // create cube potential that is spanned from (-1,-1,-1) to (1, 1, 1)
        kernel->getKernelContext().registerPotential(std::move(cubePot), t);
    }

    auto ppObs = kernel->createObservable<readdy::model::observables::ParticlePosition>(1);
    readdy::model::Vec3 lowerBound{-2.5, -2.5, -2.5}, upperBound{2.5, 2.5, 2.5};
    ppObs->setCallback([lowerBound, upperBound](readdy::model::observables::ParticlePosition::result_t currentResult) {
        readdy::model::Vec3 avg{0, 0, 0};
        bool allWithinBounds = true;
        for (auto &&v : currentResult) {
            allWithinBounds &= v >= lowerBound && v <= upperBound;
            avg += v;
        }
        avg /= currentResult.size();
        if (!allWithinBounds) readdy::log::console()->debug("Average position: {}", avg);
        ASSERT_TRUE(allWithinBounds);
    });
    auto connection = kernel->connectObservable(ppObs.get());
    run(*kernel);
}

TEST_P(TestPotentials, TestParticleStayInSphere) {
    kernel->getKernelContext().setBoxSize(10, 10, 10);
    kernel->getKernelContext().setTimeStep(.005);

    setupParticles(*kernel);

    std::array<std::string, 2> types{{"A", "B"}};
    for (auto t : types) {
        auto spherePot = kernel->createPotentialAs<readdy::model::potentials::SpherePotential>();
        spherePot->setOrigin(readdy::model::Vec3(0, 0, 0));
        spherePot->setRadius(3);
        spherePot->setForceConstant(20);
        kernel->getKernelContext().registerPotential(std::move(spherePot), t);
    }
    auto ppObs = kernel->createObservable<readdy::model::observables::ParticlePosition>(1);
    const double maxDistFromOrigin = 4.0; // at kbt=1 and force_const=20 the RMSD in a well potential would be ~0.2
    const double maxDistFromOriginSquared = maxDistFromOrigin * maxDistFromOrigin;
    ppObs->setCallback([maxDistFromOriginSquared](readdy::model::observables::ParticlePosition::result_t currentResult) {
        readdy::model::Vec3 avg{0, 0, 0};
        bool allWithinBounds = true;
        for (auto &&v : currentResult) {
            const double distanceFromOriginSquared = v * v;
            allWithinBounds &= distanceFromOriginSquared < maxDistFromOriginSquared;
            avg += v;
        }
        avg /= currentResult.size();
        if (!allWithinBounds) readdy::log::console()->debug("Average position: {}", avg);
        ASSERT_TRUE(allWithinBounds);
    });
    auto connection = kernel->connectObservable(ppObs.get());
    run(*kernel);
}

INSTANTIATE_TEST_CASE_P(TestPotentials, TestPotentials,
                        ::testing::ValuesIn(readdy::testing::getKernelsToTest()));
}
