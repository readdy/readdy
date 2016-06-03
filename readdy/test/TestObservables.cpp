/**
 * << detailed description >>
 *
 * @file TestObservables.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 02.05.16
 */

#include <boost/algorithm/string.hpp>
#include <readdy/plugin/KernelProvider.h>
#include <readdy/Simulation.h>
#include "gtest/gtest.h"

namespace m = readdy::model;

namespace {
    class  TestObservables : public ::testing::Test {
    protected:
        std::unique_ptr<m::Kernel> kernel;

        TestObservables() {
            // if we're in conda
            const char *env = std::getenv("PREFIX");
            std::string pluginDir = "lib/readdy_plugins";
            if (env) {
                auto _env = std::string(env);
                if (!boost::algorithm::ends_with(env, "/")) {
                    _env = _env.append("/");
                }
                pluginDir = _env.append(pluginDir);
            }
            readdy::plugin::KernelProvider::getInstance().loadKernelsFromDirectory(pluginDir);
            kernel = readdy::plugin::KernelProvider::getInstance().create("SingleCPU");
        }
    };

    TEST_F(TestObservables, Foo) {
        readdy::model::_internal::ObservableFactory obsf(kernel.get());
        auto x = obsf.create<m::ParticlePositionObservable>();
    }

    TEST_F(TestObservables, TestParticlePositions) {
        const unsigned int n_particles = 100;
        const double diffusionConstant = 1;
        kernel->getKernelContext().setDiffusionConstant("type", diffusionConstant);
        const double timeStep = 1.0;
        kernel->getKernelContext().setTimeStep(timeStep);
        const auto particleTypeId = kernel->getKernelContext().getParticleTypeID("type");
        const auto particles = std::vector<m::Particle>(n_particles, m::Particle(0,0,0, particleTypeId));
        kernel->getKernelStateModel().addParticles(particles);
        auto&& obs = kernel->createObservable<m::ParticlePositionObservable>();
        obs->setStride(3);
        auto &&connection = kernel->registerObservable(obs.get());

        auto&& diffuseProgram = kernel->createProgram("Diffuse");
        for(readdy::model::time_step_type t = 0; t < 100; t++) {
            diffuseProgram->execute();
            kernel->getKernelStateModel().updateModel(t, false, false);
        }

        const auto&& result = obs->getResult();

        int j = 0;
        for(auto&& p : *result) {
            BOOST_LOG_TRIVIAL(debug) << "foo " << ++j << " / " << result->size() << " (" << particles.size() << ")";
        }
        connection.disconnect();
    }

    TEST_F(TestObservables, TestCombinerObservable) {
        auto&& o1 = kernel->createObservable<m::ParticlePositionObservable>();
        auto&& o2 = kernel->createObservable<m::ParticlePositionObservable>();
        auto&& o3 = kernel->createObservable<m::TestCombinerObservable>(o1.get(), o2.get());
        auto&& connection = kernel->registerObservable(o3.get());
        auto&& diffuseProgram = kernel->createProgram("Diffuse");
        for(readdy::model::time_step_type t = 0; t < 100; t++) {
            diffuseProgram->execute();
            kernel->getKernelStateModel().updateModel(t, false, false);
        }

        const auto&& result = o3->getResult();
        for(auto&& p : *result) {
            // todo
        }

        connection.disconnect();
    }
}

