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
 * @file CPUKernelTest.cpp.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 23.06.16
 */

#include <gtest/gtest.h>
#include <readdy/plugin/KernelProvider.h>
#include <readdy/model/actions/Actions.h>
#include <readdy/model/RandomProvider.h>

namespace {

TEST(CPUTestKernel, TestKernelLoad) {
    auto kernel = readdy::plugin::KernelProvider::getInstance().create("CPU");

    kernel->getKernelContext().setBoxSize(10, 10, 10);
    kernel->getKernelContext().registerParticleType("X", .55, 1.);
    kernel->getKernelContext().setPeriodicBoundary(true, true, true);
    kernel->registerReaction<readdy::model::reactions::Decay>("X decay", "X", .5);
    kernel->registerReaction<readdy::model::reactions::Fission>("X fission", "X", "X", "X", .00, .5);

    auto &&integrator = kernel->createAction<readdy::model::actions::EulerBDIntegrator>(1);
    auto &&neighborList = kernel->createAction<readdy::model::actions::UpdateNeighborList>();
    auto &&forces = kernel->createAction<readdy::model::actions::CalculateForces>();
    auto &&reactions = kernel->createAction<readdy::model::actions::reactions::GillespieParallel>(1);

    auto pp_obs = kernel->createObservable<readdy::model::observables::Positions>(1);
    auto connection = kernel->connectObservable(pp_obs.get());

    const int n_particles = 500;
    const auto typeId = kernel->getKernelContext().getParticleTypeID("X");
    std::vector<readdy::model::Particle> particlesToBeginWith{n_particles, {0, 0, 0, typeId}};
    readdy::log::console()->debug("n_particles={}", particlesToBeginWith.size());
    kernel->getKernelStateModel().addParticles(particlesToBeginWith);

    kernel->getKernelContext().configure();

    neighborList->perform();
    for (size_t t = 0; t < 20; t++) {

        forces->perform();
        integrator->perform();
        neighborList->perform();

        reactions->perform();

        pp_obs->evaluate();
        readdy::log::console()->debug("\tcurrently n particles: {}", pp_obs->getResult().size());
    }

    connection.disconnect();
}
}