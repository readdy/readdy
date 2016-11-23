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
#include <readdy/common/Timer.h>
#include <readdy/model/Utils.h>

namespace {

void runPerformanceTest(readdy::model::Kernel &kernel, readdy::model::observables::time_step_type steps = 3) {
    using timer = readdy::util::Timer;

    auto stdRand = [](double lower = 0.0, double upper = 1.0) -> double {
        return static_cast <double> (std::rand()) / (RAND_MAX / (upper - lower)) + lower;
    };

    kernel.getKernelContext().setBoxSize(30.0, 80.0, 30.0);
    kernel.getKernelContext().setKBT(1.0);
    kernel.getKernelContext().setTimeStep(.001);
    kernel.getKernelContext().setPeriodicBoundary(false, false, false);

    std::string types[]{"A", "B", "C"};
    std::string tup[3][2] = {{"A", "B"}, {"A", "C"}, {"B", "C"}};
    {
        double x = 1.0;
        for (auto &&t : types) {
            kernel.getKernelContext().setDiffusionConstant(t, x++);
            kernel.getKernelContext().setParticleRadius(t, 1.0);
        }
    }


    for(auto i = 0; i < 3; ++i){
        auto t = tup[i];
        auto repulsion = kernel.createPotentialAs<readdy::model::potentials::HarmonicRepulsion>();
        repulsion->setForceConstant(1.0);
        kernel.getKernelContext().registerPotential(std::move(repulsion), t[0], t[1]);
    }
    {
        const auto verts = kernel.getKernelContext().getBoxBoundingVertices();
        int i = 0;
        for(auto type : types) {
            auto box = kernel.createPotentialAs<readdy::model::potentials::CubePotential>();
            box->setForceConstant(1.0);
            box->setOrigin(std::get<0>(verts) + readdy::model::Vec3(1, 1, 1));
            box->setExtent((std::get<1>(verts) - std::get<0>(verts)) - 2);
            kernel.getKernelContext().registerPotential(std::move(box), type);
        }
    }

    const unsigned int nParticles = 50000;
    for (unsigned long _ = 0; _ < nParticles; ++_) {
        for (const auto &t : types) {
            readdy::model::Particle p{stdRand(-15, 15), stdRand(-40, 40), stdRand(-15, 15),
                                      kernel.getKernelContext().getParticleTypeID(t)};
            kernel.getKernelStateModel().addParticle(p);
        }
    }

    kernel.registerReaction<readdy::model::reactions::Fusion>("A+B->C", "A", "B", "C", .05, .05);
    kernel.registerReaction<readdy::model::reactions::Fission>("C->A+B", "C", "A", "B", .05, .05);

    auto &&integrator = kernel.createProgram<readdy::model::programs::EulerBDIntegrator>();
    auto &&neighborList = kernel.createProgram<readdy::model::programs::UpdateNeighborList>();
    neighborList->setSkinSize(20*readdy::model::util::getMaximumDisplacement(kernel.getKernelContext()));
    auto &&forces = kernel.createProgram<readdy::model::programs::CalculateForces>();
    auto &&reactionsProgram = kernel.createProgram<readdy::model::programs::reactions::GillespieParallel>();

    auto obs = kernel.createObservable<readdy::model::NParticlesObservable>(0);
    obs->setCallback([](const std::vector<unsigned long> n) {
        readdy::log::console()->debug("have n particles = {}", n[0]);
    });
    auto connection = kernel.connectObservable(obs.get());

    double t_forces = 0, t_integrator = 0, t_nl = 0, t_reactions = 0;

    kernel.getKernelContext().configure();

    {
        timer c("neighbor list init");
        neighborList->execute();
        t_nl += c.getSeconds();
    }
    for (readdy::model::observables::time_step_type t = 0; t < steps; ++t) {
        readdy::log::console()->debug("----------");
        readdy::log::console()->debug("t = {}", t);
        kernel.evaluateObservables(t);
        {
            timer c("forces");
            forces->execute();
            t_forces += c.getSeconds();
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

        {
            timer c("neighbor list");
            neighborList->execute();
            t_nl += c.getSeconds();
        }
    }

    std::cout << "--------------------------------------------------------------" << std::endl;
    std::cout << "Average time for calculating forces: " << t_forces / steps << std::endl;
    std::cout << "Average time for the integrator:     " << t_integrator / steps << std::endl;
    std::cout << "Average time for the neighbor list:  " << t_nl / (2*steps+1) << std::endl;;
    std::cout << "Average time for handling reactions: " << t_reactions / steps << std::endl;;
    std::cout << "--------------------------------------------------------------" << std::endl;;
}

TEST(TestPerformance, SingleCPU) {
    //auto kernel = readdy::plugin::KernelProvider::getInstance().create("SingleCPU");
    //runPerformanceTest(*kernel);
}

TEST(TestPerformance, CPU) {
    auto kernel = readdy::plugin::KernelProvider::getInstance().create("CPU");
    runPerformanceTest(*kernel, 20);
}

}