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


#include "gtest/gtest.h"
#include <readdy/api/Simulation.h>

using namespace readdy;

namespace {

TEST(TestSimulation, TestObservables) {
    Simulation simulation{"SingleCPU"};
    simulation.context().boxSize() = {{10, 10, 10}};
    unsigned int n_particles = 103;
    readdy::scalar diffusionConstant = 1;
    simulation.context().particleTypes().add("type", diffusionConstant);
    for (auto _ = 0; _ < n_particles; ++_) {
        simulation.addParticle("type", 0, 0, 0);
    }
    readdy::scalar timestep = 1;

    int n_callbacks = 0;
    simulation.registerObservable(simulation.observe().positions(1),
                                  [&n_callbacks](const readdy::model::observables::Positions::result_type &result) {
                                      ++n_callbacks;
                                      EXPECT_EQ(103, result.size());
                                  });
    simulation.run(100, timestep);
    EXPECT_EQ(101, n_callbacks);
}
}
