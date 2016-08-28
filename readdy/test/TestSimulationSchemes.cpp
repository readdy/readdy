/**
 * << detailed description >>
 *
 * @file TestSimulationSchemes.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 23.08.16
 */

#include <gtest/gtest.h>
#include <readdy/plugin/KernelProvider.h>
#include <readdy/Simulation.h>

namespace api = readdy::api;

namespace {

    TEST(TestSchemes, Sanity) {
        auto kernel = readdy::plugin::KernelProvider::getInstance().create("SingleCPU");
        api::SchemeConfigurator<api::ReaDDyScheme> c(kernel.get());
        c.evaluateObservables(false)
                .includeForces(false)
                .withIntegrator(kernel->createProgram<readdy::model::programs::EulerBDIntegrator>())
                .configure()->run(10);
    }

    TEST(TestSchemes, SimulationObject) {
        readdy::Simulation sim;
        sim.setKernel("SingleCPU");
        sim.setTimeStep(.5);
        sim.setBoxSize(1, 1, 1);
        sim.runScheme().configureAndRun(5);

        /**
         * use ReaDDyScheme without defaults
         */
        sim.runScheme<readdy::api::ReaDDyScheme>(false)
                .withIntegrator<readdy::model::programs::EulerBDIntegrator>()
                .withReactionScheduler<readdy::model::programs::reactions::UncontrolledApproximation>()
                .configureAndRun(100);

        /**
         * default: readdy scheme, use defaults = true
         */
        sim.runScheme()
                .withIntegrator<readdy::model::programs::EulerBDIntegrator>()
                .withReactionScheduler<readdy::model::programs::reactions::UncontrolledApproximation>()
                .configureAndRun(100);
    }

}