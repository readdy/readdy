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
 * @date 23.08.16
 */

#include <gtest/gtest.h>
#include <readdy/plugin/KernelProvider.h>
#include <readdy/Simulation.h>

namespace api = readdy::api;

namespace {

TEST(TestSchemes, Sanity) {
    auto kernel = readdy::plugin::KernelProvider::getInstance().create("SingleCPU");
    api::SchemeConfigurator <api::ReaDDyScheme> c(kernel.get());
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