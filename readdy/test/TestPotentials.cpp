/**
 * << detailed description >>
 *
 * @file TestPotentials.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 27.06.16
 */

#include <gtest/gtest.h>
#include <readdy/common/Utils.h>
#include <readdy/plugin/KernelProvider.h>
#include <readdy/model/potentials/PotentialsOrder2.h>
#include <readdy/model/potentials/PotentialsOrder1.h>
#include <readdy/model/programs/Programs.h>

namespace {

    TEST(TestPotentials, TestParticlesStayInBox) {
        auto kernel = readdy::plugin::KernelProvider::getInstance().create("SingleCPU");
        kernel->getKernelContext().setDiffusionConstant("A", 1);
        kernel->getKernelContext().setDiffusionConstant("B", .1);
        kernel->getKernelContext().setParticleRadius("A", .1);
        kernel->getKernelContext().setParticleRadius("B", .01);
        kernel->getKernelContext().setPeriodicBoundary(false, false, false);
        kernel->getKernelContext().setBoxSize(5, 5, 5);
        kernel->getKernelContext().setTimeStep(.005);
        const unsigned int nParticlesA = 10;
        const unsigned int nParticlesB = 10;

        for (auto i = 0; i < nParticlesA; ++i) kernel->addParticle("A", {0, 0, 0});
        for (auto i = 0; i < nParticlesB; ++i) kernel->addParticle("B", {0, 0, 0});

        auto harmonicRepulsion = kernel->createPotentialAs<readdy::model::potentials::HarmonicRepulsion>();
        harmonicRepulsion->setForceConstant(.01);
        kernel->getKernelContext().registerOrder2Potential(harmonicRepulsion.get(), "A", "B");

        // create cube potential that is spanned from (-1,-1,-1) to (1, 1, 1)
        auto cubePot = kernel->createPotentialAs<readdy::model::potentials::CubePotential>();
        cubePot->setOrigin({-1, -1, -1});
        cubePot->setConsiderParticleRadius(true);
        cubePot->setExtent({2, 2, 2});
        cubePot->setForceConstant(10);

        kernel->getKernelContext().registerOrder1Potential(cubePot.get(), "A");
        kernel->getKernelContext().registerOrder1Potential(cubePot.get(), "B");
        auto ppObs = kernel->createObservable<readdy::model::ParticlePositionObservable>(1);
        readdy::model::Vec3 lowerBound{-2.5, -2.5, -2.5}, upperBound{2.5, 2.5, 2.5};
        ppObs->setCallback([lowerBound, upperBound](readdy::model::ParticlePositionObservable::result_t currentResult) {
            readdy::model::Vec3 avg{0, 0, 0};
            bool allWithinBounds = true;
            for (auto &&v : currentResult) {
                allWithinBounds &= v >= lowerBound && v <= upperBound;
                avg += v;
            }
            avg /= currentResult.size();
            if (!allWithinBounds) BOOST_LOG_TRIVIAL(debug) << "Average position: " << avg;
            ASSERT_TRUE(allWithinBounds);
        });
        auto connection = kernel->connectObservable(ppObs.get());

        unsigned int nSteps = 1000;
        auto &&integrator = kernel->createProgram<readdy::model::programs::EulerBDIntegrator>();
        auto &&nl = kernel->createProgram<readdy::model::programs::UpdateNeighborList>();
        auto &&forces = kernel->createProgram<readdy::model::programs::CalculateForces>();
        nl->execute();
        for (readdy::model::time_step_type &&t = 0; t < nSteps; ++t) {
            forces->execute();
            integrator->execute();
            nl->execute();
        }
    }
}