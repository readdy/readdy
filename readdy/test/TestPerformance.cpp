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
#include <readdy/testing/Timer.h>

namespace {

    void runPerformanceTest(readdy::model::Kernel &kernel, readdy::model::time_step_type steps = 3) {
        using timer = readdy::testing::Timer;

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
        auto &&neighborList = kernel.createProgram<readdy::model::programs::UpdateNeighborList>();
        auto &&forces = kernel.createProgram<readdy::model::programs::CalculateForces>();
        auto &&reactionsProgram = kernel.createProgram<readdy::model::programs::reactions::Gillespie>();
        kernel.getKernelContext().configure();

        auto obs = kernel.createObservable<readdy::model::NParticlesObservable>(0);
        obs->setCallback([] (const std::vector<unsigned long> n) {
            BOOST_LOG_TRIVIAL(debug) << "have n particles = " << n[0];
        });
        auto connection = kernel.connectObservable(obs.get());

        double t_forces=0, t_integrator=0, t_nl=0, t_reactions = 0;

        neighborList->execute();
        for(readdy::model::time_step_type t = 0; t < steps; ++t) {
            BOOST_LOG_TRIVIAL(debug) << "----------";
            BOOST_LOG_TRIVIAL(debug) << "t = " << t;
            kernel.evaluateObservables(t);
            {
                timer c("forces");
                forces->execute();
                t_forces+=c.getSeconds();
            }
            {
                timer c("integrator");
                integrator->execute();
                t_integrator += c.getSeconds();
            }
            {
                timer c("neighbor list");
                neighborList->execute();
                t_nl += c.getSeconds();
            }

            {
                timer c("reactions");
                reactionsProgram->execute();
                t_reactions += c.getSeconds();
            }
        }

        BOOST_LOG_TRIVIAL(debug) << "--------------------------------------------------------------";
        BOOST_LOG_TRIVIAL(debug) << "Average time for calculating forces: " << t_forces / steps;
        BOOST_LOG_TRIVIAL(debug) << "Average time for the integrator:     " << t_integrator / steps;
        BOOST_LOG_TRIVIAL(debug) << "Average time for the neighbor list:  " << t_nl / steps;
        BOOST_LOG_TRIVIAL(debug) << "Average time for handling reactions: " << t_reactions / steps;
        BOOST_LOG_TRIVIAL(debug) << "--------------------------------------------------------------";
    }

    TEST(TestPerformance, SingleCPU) {
        {
            auto kernel = readdy::plugin::KernelProvider::getInstance().create("SingleCPU");
            runPerformanceTest(*kernel);
        }
    }

    TEST(TestPerformance, CPU) {
        {
            auto kernel = readdy::plugin::KernelProvider::getInstance().create("CPU");
            runPerformanceTest(*kernel);
        }
    }

}