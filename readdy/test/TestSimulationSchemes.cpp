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

    explicit TestSchemes() : simulation(GetParam()) {}
};

TEST_P(TestSchemes, CorrectNumberOfTimesteps) {
    unsigned int counter = 0;
    auto increment = [&counter](readdy::model::observables::NParticles::result_type result) {
        counter++;
    };
    auto obsHandle = simulation.registerObservable(simulation.observe().nParticles(1), increment);
    simulation.run(3, 0.1);
    EXPECT_EQ(counter, 4);
}

TEST_P(TestSchemes, StoppingCriterionSimple) {
    unsigned int counter = 0;
    auto increment = [&counter](readdy::model::observables::NParticles::result_type result) {
        counter++;
    };
    auto obsHandle = simulation.registerObservable(simulation.observe().nParticles(1), increment);
    auto shallContinue = [](readdy::time_step_type currentStep) {
        return currentStep < 5;
    };
    auto loop = simulation.createLoop(.1);
    loop.run(shallContinue);
    EXPECT_EQ(counter, 6);
}

TEST_P(TestSchemes, ComplexStoppingCriterion) {
    simulation.context().particleTypes().add("A", 0.);
    // A -> A + A, with probability = 1 each timestep. After 3 timesteps there will be 8 particles. The counter will be 4 by then.
    simulation.context().reactions().addFission("bla", "A", "A", "A", 1e8, 0.);
    simulation.addParticle("A", 0, 0, 0);
    unsigned int counter = 0;
    bool doStop = false;
    auto increment = [&counter, &doStop](readdy::model::observables::NParticles::result_type result) {
        counter++;
        if (result[0] >= 8) {
            doStop = true;
        }
    };
    auto obsHandle = simulation.registerObservable(simulation.observe().nParticles(1), increment);
    auto shallContinue = [&doStop](readdy::time_step_type currentStep) {
        return !doStop;
    };
    simulation.createLoop(1.).run(shallContinue);
    EXPECT_EQ(counter, 4);
}

TEST_P(TestSchemes, SkinSizeSanity) {
    simulation.context().particleTypes().add("A", 1.);
    simulation.context().boxSize() = {{10., 10., 10.}};
    simulation.context().periodicBoundaryConditions() = {{true, true, true}};
    simulation.context().potentials().addHarmonicRepulsion("A", "A", 1., 2.);
    simulation.addParticle("A", 0., 0., 0.);
    simulation.addParticle("A", 1.5, 0., 0.);
    auto loop = simulation.createLoop(.001);
    loop.skinSize() = 1.;
    loop.run(10);
}

INSTANTIATE_TEST_CASE_P(TestSchemesCore, TestSchemes, ::testing::ValuesIn(readdy::testing::getKernelsToTest()));

}