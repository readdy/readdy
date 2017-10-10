/********************************************************************
 * Copyright © 2016 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * This file is part of ReaDDy.                                     *
 *                                                                  *
 * ReaDDy is free software: you can redistribute it and/or modify   *
 * it under the terms of the GNU Lesser General Public License as   *
 * published by the Free Software Foundation, either version 3 of   *
 * the License, or (at your option) any later version.              *
 *                                                                  *
 * This program is distributed in the hope that it will be useful,  *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of   *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    *
 * GNU Lesser General Public License for more details.              *
 *                                                                  *
 * You should have received a copy of the GNU Lesser General        *
 * Public License along with this program. If not, see              *
 * <http://www.gnu.org/licenses/>.                                  *
 ********************************************************************/


/**
 * << detailed description >>
 *
 * @file TestSimulationSchemes.cpp
 * @brief << brief description >>
 * @author clonker
 * @author chrisfro**
 * @date 23.08.16
 */

#include <gtest/gtest.h>
#include <readdy/plugin/KernelProvider.h>
#include <readdy/api/Simulation.h>
#include <readdy/testing/Utils.h>

namespace api = readdy::api;

namespace {

class TestSchemes : public ::testing::TestWithParam<std::string> {
public:
    readdy::Simulation simulation;

    explicit TestSchemes() {
        simulation.setKernel(GetParam());
    }
};

TEST_P(TestSchemes, SimulationObject) {
    simulation.setBoxSize(1, 1, 1);
    simulation.runScheme().configureAndRun(5, .5);

    /**
     * use ReaDDyScheme without defaults
     */
    simulation.runScheme<readdy::api::ReaDDyScheme>(false)
            .withIntegrator<readdy::model::actions::EulerBDIntegrator>()
            .withReactionScheduler<readdy::model::actions::reactions::UncontrolledApproximation>()
            .configureAndRun(5, .5);

    /**
     * default: readdy scheme, use defaults = true
     */
    simulation.runScheme()
            .includeForces(false)
            .withIntegrator<readdy::model::actions::EulerBDIntegrator>()
            .withReactionScheduler<readdy::model::actions::reactions::UncontrolledApproximation>()
            .configureAndRun(5, .5);

    /**
     * use AdvancedScheme
     */
    simulation.runScheme<readdy::api::AdvancedScheme>(false)
            .includeForces(true)
            .withIntegrator<readdy::model::actions::EulerBDIntegrator>()
            .includeCompartments(true)
            .configureAndRun(5, .5);
}

TEST_P(TestSchemes, CorrectNumberOfTimesteps) {
    unsigned int counter = 0;
    auto increment = [&counter](readdy::model::observables::NParticles::result_type result) {
        counter++;
    };
    auto obsHandle = simulation.registerObservable<readdy::model::observables::NParticles>(increment, 1);
    simulation.run(3, 0.1);
    EXPECT_EQ(counter, 4);
}

TEST_P(TestSchemes, StoppingCriterionSimple) {
    unsigned int counter = 0;
    auto increment = [&counter](readdy::model::observables::NParticles::result_type result) {
        counter++;
    };
    auto obsHandle = simulation.registerObservable<readdy::model::observables::NParticles>(increment, 1);
    auto shallContinue = [](readdy::time_step_type currentStep) {
        return currentStep < 5;
    };
    auto scheme = simulation.runScheme(true).configure(0.1);
    scheme->run(shallContinue);
    EXPECT_EQ(counter, 6);
}

TEST_P(TestSchemes, ComplexStoppingCriterion) {
    simulation.registerParticleType("A", 0.);
    // A -> A + A, with probability = 1 each timestep. After 3 timesteps there will be 8 particles. The counter will be 4 by then.
    simulation.registerFissionReaction("bla", "A", "A", "A", 1000., 0.);
    simulation.addParticle("A", 0, 0, 0);
    unsigned int counter = 0;
    bool doStop = false;
    auto increment = [&counter, &doStop](readdy::model::observables::NParticles::result_type result) {
        counter++;
        if (result[0] >= 8) {
            doStop = true;
        }
    };
    auto obsHandle = simulation.registerObservable<readdy::model::observables::NParticles>(increment, 1);
    auto shallContinue = [&doStop](readdy::time_step_type currentStep) {
        return !doStop;
    };
    auto scheme = simulation.runScheme(true).configure(1.);
    scheme->run(shallContinue);
    EXPECT_EQ(counter, 4);
}

TEST_P(TestSchemes, SkinSizeSanity) {
    simulation.registerParticleType("A", 1.);
    simulation.setBoxSize(10., 10., 10.);
    simulation.setPeriodicBoundary({true, true, true});
    simulation.registerHarmonicRepulsionPotential("A", "A", 1., 2.);
    simulation.addParticle("A", 0., 0., 0.);
    simulation.addParticle("A", 1.5, 0., 0.);
    readdy::api::SchemeConfigurator<readdy::api::ReaDDyScheme> configurator = simulation.runScheme(true);
    configurator.withSkinSize(1.);
    auto scheme = configurator.configure(0.001);
    scheme->run(10);
}

INSTANTIATE_TEST_CASE_P(TestSchemesCore, TestSchemes, ::testing::ValuesIn(readdy::testing::getKernelsToTest()));

}