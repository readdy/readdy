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
        kernel->getKernelContext().registerDeathReaction("X decay", "X", .05);
        kernel->getKernelContext().registerFissionReaction("X fission", "X", "X", "X", .5, .00);

        auto &&diffuseProgram = kernel->createProgram<readdy::model::programs::DiffuseProgram>();
        auto &&updateModelProgram = kernel->createProgram<readdy::model::programs::UpdateStateModelProgram>();
        auto &&reactionsProgram = kernel->createProgram<readdy::model::programs::DefaultReactionProgram>();

        auto pp_obs = kernel->createObservable<readdy::model::ParticlePositionObservable>(1);
        auto connection = kernel->connectObservable(pp_obs.get());

        const int n_particles = 2000;
        const unsigned int typeId = kernel->getKernelContext().getParticleTypeID("X");
        std::vector<readdy::model::Particle> particlesToBeginWith {n_particles, {0,0,0,typeId}};
        BOOST_LOG_TRIVIAL(debug) << "n_particles="<<particlesToBeginWith.size();
        kernel->getKernelStateModel().addParticles(particlesToBeginWith);

        std::unique_ptr<readdy::model::RandomProvider> rand = std::make_unique<readdy::model::RandomProvider>();

        for(size_t t = 0; t < 1000; t++) {

            diffuseProgram->execute();
            updateModelProgram->configure(t, false);
            updateModelProgram->execute();

            reactionsProgram->execute();
            updateModelProgram->configure(t, false);
            updateModelProgram->execute();

            pp_obs->evaluate();
            BOOST_LOG_TRIVIAL(debug) << "\tcurrently n particles: " << pp_obs->getResult().size();

        }

        connection.disconnect();
    }
}