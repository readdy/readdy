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
 * @file TestAggregators.cpp
 * @brief Test correctness of aggregators
 * @author chrisfroe
 * @date 07.11.16
 */

#include <gmock/gmock.h>
#include <readdy/testing/KernelTest.h>
#include <readdy/testing/Utils.h>

namespace {

class TestAggregators : public KernelTest {
};

TEST_P(TestAggregators, TestMeanSquaredDisplacement) {
    kernel->context().particleTypes().add("A", 1.);
    for (auto i=0; i<5; ++i) kernel->addParticle("A", readdy::Vec3(0, 0, 0));
    auto obs = kernel->observe().particles(1);
    auto msd = kernel->observe().msd(1, std::vector<std::string>(), obs.get());
    auto connection = kernel->connectObservable(msd.get());
    for (auto i=0; i<3; ++i) {
        kernel->evaluateObservables(0);
        const auto& result = msd->getResult();
        const auto& resultMsd = std::get<1>(result);
        const auto& resultTime = std::get<0>(result);
        EXPECT_EQ(resultMsd.size(), 1) << "repetitive evaluation at t=0 should always yield result of length 1";
        EXPECT_EQ(resultTime.size(), 1);
        EXPECT_THAT(resultMsd, ::testing::ElementsAre(0)) << "particles have not moved -> msd=0";
        EXPECT_THAT(resultTime, ::testing::ElementsAre(0));
    }
    for (auto i=0; i<3; ++i) {
        kernel->evaluateObservables(1);
        const auto& result = msd->getResult();
        const auto& resultMsd = std::get<1>(result);
        const auto& resultTime = std::get<0>(result);
        EXPECT_EQ(resultMsd.size(), 2) << "repetitive evaluation at t=1 should always yield result of length 2";
        EXPECT_EQ(resultTime.size(), 2);
        EXPECT_THAT(resultMsd, ::testing::ElementsAre(0, 0)) << "particles have not moved -> msd=0";
        EXPECT_THAT(resultTime, ::testing::ElementsAre(0, 1));
    }
}

TEST_P(TestAggregators, TestTrivial) {
    kernel->context().particleTypes().add("A", 1.);
    for (auto i=0; i<5; ++i) kernel->addParticle("A", readdy::Vec3(4, 2, 0));
    auto obs = kernel->observe().positions(1);
    auto traj = kernel->observe().collect(1, obs.get());
    auto connection = kernel->connectObservable(traj.get());
    for (auto i=0; i<3; ++i) {
        kernel->evaluateObservables(0);
        const auto& result = traj->getResult();
        const auto& resultPos = std::get<1>(result);
        const auto& resultTime = std::get<0>(result);
        EXPECT_EQ(resultTime.size(), 1);
        EXPECT_EQ(resultPos.size(), 1);
        EXPECT_THAT(resultTime, ::testing::ElementsAre(0));
        EXPECT_EQ(resultPos[0].size(), 5);
        for (auto j=0; j<5; ++j) {
            EXPECT_VEC3_EQ(resultPos[0][j], readdy::Vec3(4, 2, 0));
        }
    }
    for (auto i=0; i<3; ++i) {
        kernel->evaluateObservables(1);
        const auto& result = traj->getResult();
        const auto& resultPos = std::get<1>(result);
        const auto& resultTime = std::get<0>(result);
        EXPECT_EQ(resultTime.size(), 2);
        EXPECT_EQ(resultPos.size(), 2);
        EXPECT_THAT(resultTime, ::testing::ElementsAre(0, 1));
        EXPECT_EQ(resultPos[1].size(), 5);
        for (auto j=0; j<5; ++j) {
            EXPECT_VEC3_EQ(resultPos[1][j], readdy::Vec3(4, 2, 0));
        }
    }
}

INSTANTIATE_TEST_CASE_P(TestAggregators, TestAggregators, ::testing::ValuesIn(readdy::testing::getKernelsToTest()));

}
