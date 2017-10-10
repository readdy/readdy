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

struct MSDAggregator {

    readdy::scalar msd = 0;
    long T = 0;
    std::vector<readdy::Vec3> initialPositions;
    std::shared_ptr<readdy::scalar> result = std::make_shared<readdy::scalar>(0);

    void operator()(readdy::model::observables::Positions::result_type positions) {
        auto it_init = initialPositions.begin();
        auto it_pos = positions.begin();
        while (it_pos != positions.end()) {
            msd += (*it_init - *it_pos) * (*it_init - *it_pos);
            ++it_init;
            ++it_pos;
        }
        ++T;
        *result = msd / T;
    }
};

class TestSimulation : public ::testing::Test {
protected:
    Simulation simulation;

    TestSimulation() {
        simulation.setKernel("SingleCPU");
    }
};

TEST_F(TestSimulation, TestKBT) {
    simulation.setKBT(42);
    EXPECT_EQ(42, simulation.getKBT());
}

TEST_F(TestSimulation, TestPeriodicBdry) {
    simulation.setPeriodicBoundary({true, false, true});
    auto boundary = simulation.getPeriodicBoundary();
    EXPECT_TRUE(boundary[0]);
    EXPECT_FALSE(boundary[1]);
    EXPECT_TRUE(boundary[2]);
}

TEST_F(TestSimulation, TestBoxSize) {
    simulation.setBoxSize(10, 11, 12);
    auto box_size = simulation.getBoxSize();
    EXPECT_EQ(box_size[0], 10);
    EXPECT_EQ(box_size[1], 11);
    EXPECT_EQ(box_size[2], 12);
}

TEST_F(TestSimulation, TestMeanSquaredDisplacement) {
    Simulation simulation;
    simulation.setKernel("SingleCPU");
    simulation.setBoxSize(10, 10, 10);
    unsigned int n_particles = 103;
    readdy::scalar diffusionConstant = 1;
    simulation.registerParticleType("type", diffusionConstant);
    for (auto _ = 0; _ < n_particles; ++_) {
        simulation.addParticle("type", 0, 0, 0);
    }
    readdy::scalar timestep = 1;
    MSDAggregator aggregator;
    aggregator.initialPositions = simulation.getAllParticlePositions();
    simulation.registerObservable<readdy::model::observables::Positions>(aggregator, 1);
    simulation.run(100, timestep);
    auto positions = simulation.getAllParticlePositions();
    readdy::scalar msd = 0;
    for (auto &&position : positions) {
        msd += position * position;
    }
    msd /= positions.size();
    readdy::log::debug("mean squared displacement: {}", msd);
    readdy::log::debug("mean squared displacement2: {}", *aggregator.result);
}

TEST_F(TestSimulation, TestObservables) {
    Simulation simulation;
    simulation.setKernel("SingleCPU");
    simulation.setBoxSize(10, 10, 10);
    unsigned int n_particles = 103;
    readdy::scalar diffusionConstant = 1;
    simulation.registerParticleType("type", diffusionConstant);
    for (auto _ = 0; _ < n_particles; ++_) {
        simulation.addParticle("type", 0, 0, 0);
    }
    readdy::scalar timestep = 1;

    int n_callbacks = 0;
    simulation.registerObservable<readdy::model::observables::Positions>(
            [&n_callbacks](const readdy::model::observables::Positions::result_type &result) -> void {
                ++n_callbacks;
                EXPECT_EQ(103, result.size());
            }, 1);
    simulation.run(100, timestep);
    EXPECT_EQ(101, n_callbacks);
}
}