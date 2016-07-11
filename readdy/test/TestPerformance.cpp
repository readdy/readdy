/**
 * This file includes a larger simulation with fluctuating particle numbers so that the performance can be
 * evaluated and bottlenecks be identified.
 *
 * @file TestPerformance.cpp
 * @brief Test including a larger simulation with varying particle numbers. Per default not included in tests.
 * @author clonker
 * @date 11.07.16
 */

#include <gtest/gtest.h>
#include <readdy/plugin/KernelProvider.h>
#include <readdy/model/programs/Programs.h>
#include <readdy/model/potentials/PotentialsOrder2.h>

namespace {

    void runPerformanceTest(readdy::model::Kernel &kernel, readdy::model::time_step_type steps = 100) {
        std::srand((unsigned int) std::time(0));

        auto stdRand = [] (double lower = 0.0, double upper = 1.0) -> double {
            return static_cast <double> (std::rand()) / (RAND_MAX / (upper - lower)) + lower;
        };

        kernel.getKernelContext().setBoxSize(15.0, 15.0, 15.0);
        kernel.getKernelContext().setKBT(1.0);
        kernel.getKernelContext().setTimeStep(1.);

        std::string types [] {"A", "B", "C"};
        {
            double x = 1.0;
            for (auto &&t : types) {
                kernel.getKernelContext().setDiffusionConstant(t, x++);
                kernel.getKernelContext().setParticleRadius(t, 1.0);
            }
        }

        {
            auto repulsion = kernel.createPotentialAs<readdy::model::potentials::HarmonicRepulsion>();
            repulsion->setForceConstant(1.0);
            kernel.getKernelContext().registerOrder2Potential(repulsion.get(), "A", "B");
            kernel.getKernelContext().registerOrder2Potential(repulsion.get(), "A", "C");
            kernel.getKernelContext().registerOrder2Potential(repulsion.get(), "B", "C");
        }

        const unsigned int nParticles = 1000;
        for(unsigned long _ = 0; _ < nParticles; ++_) {
            for(const auto& t : types) {
                readdy::model::Particle p{stdRand(-7.5, 7.5), stdRand(-7.5, 7.5), stdRand(-7.5, 7.5),
                                          kernel.getKernelContext().getParticleTypeID(t)};
                kernel.getKernelStateModel().addParticle(p);
            }
        }

        kernel.getKernelContext().registerFusionReaction("A+B->C", "A", "B", "C", .5, 2.0);
        kernel.getKernelContext().registerFissionReaction("C->A+B", "C", "A", "B", 2.0, .5);

        auto &&integrator = kernel.createProgram<readdy::model::programs::EulerBDIntegrator>();
        auto &&updateModelProgram = kernel.createProgram<readdy::model::programs::UpdateStateModelProgram>();
        auto &&reactionsProgram = kernel.createProgram<readdy::model::programs::DefaultReactionProgram>();
        kernel.getKernelContext().configure();

        auto obs = kernel.createObservable<readdy::model::NParticlesObservable>(10);
        obs->setCallback([] (const long n) {
            BOOST_LOG_TRIVIAL(debug) << "have n particles = " << n;
        });
        auto connection = kernel.connectObservable(obs.get());

        for(readdy::model::time_step_type t = 0; t < steps; ++t) {
            updateModelProgram->configure(t, true);
            updateModelProgram->execute();
            integrator->execute();

            updateModelProgram->configure(t, false);
            updateModelProgram->execute();
            reactionsProgram->execute();
        }
    }

    TEST(TestPerformance, SingleCPU) {
        auto kernel = readdy::plugin::KernelProvider::getInstance().create("SingleCPU");
        runPerformanceTest(*kernel);
    }

    TEST(TestPerformance, CPU) {
        auto kernel = readdy::plugin::KernelProvider::getInstance().create("CPU");
        runPerformanceTest(*kernel);
    }

}