/**
 * << detailed description >>
 *
 * @file CPUKernelTest.cpp.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 23.06.16
 */

#include <gtest/gtest.h>
#include <boost/algorithm/string.hpp>
#include <readdy/plugin/KernelProvider.h>
#include <readdy/model/programs/Programs.h>
#include <readdy/model/RandomProvider.h>

namespace {

    TEST(CPUTestKernel, TestKernelLoad) {
        auto kernel = readdy::plugin::KernelProvider::getInstance().create("CPU");

        kernel->getKernelContext().setBoxSize(10, 10, 10);
        kernel->getKernelContext().setTimeStep(1);
        kernel->getKernelContext().setDiffusionConstant("X", .55);
        kernel->registerReaction<readdy::model::reactions::Decay>("X decay", "X", .5);
        kernel->registerReaction<readdy::model::reactions::Fission>("X fission", "X", "X", "X", .00, .5);

        auto &&integrator = kernel->createProgram<readdy::model::programs::EulerBDIntegrator>();
        auto &&neighborList = kernel->createProgram<readdy::model::programs::UpdateNeighborList>();
        auto &&forces = kernel->createProgram<readdy::model::programs::CalculateForces>();
        auto &&reactionsProgram = kernel->createProgram<readdy::model::programs::reactions::UncontrolledApproximation>();

        auto pp_obs = kernel->createObservable<readdy::model::ParticlePositionObservable>(1);
        auto connection = kernel->connectObservable(pp_obs.get());

        const int n_particles = 2000;
        const unsigned int typeId = kernel->getKernelContext().getParticleTypeID("X");
        std::vector<readdy::model::Particle> particlesToBeginWith {n_particles, {0,0,0,typeId}};
        BOOST_LOG_TRIVIAL(debug) << "n_particles="<<particlesToBeginWith.size();
        kernel->getKernelStateModel().addParticles(particlesToBeginWith);

        kernel->getKernelContext().configure();

        neighborList->execute();
        for(size_t t = 0; t < 1000; t++) {

            forces->execute();
            integrator->execute();
            neighborList->execute();

            reactionsProgram->execute();

            pp_obs->evaluate();
            BOOST_LOG_TRIVIAL(debug) << "\tcurrently n particles: " << pp_obs->getResult().size();

        }

        connection.disconnect();
    }
}