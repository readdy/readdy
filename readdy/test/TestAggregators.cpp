/********************************************************************
 * Copyright © 2018 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * Redistribution and use in source and binary forms, with or       *
 * without modification, are permitted provided that the            *
 * following conditions are met:                                    *
 *  1. Redistributions of source code must retain the above         *
 *     copyright notice, this list of conditions and the            *
 *     following disclaimer.                                        *
 *  2. Redistributions in binary form must reproduce the above      *
 *     copyright notice, this list of conditions and the following  *
 *     disclaimer in the documentation and/or other materials       *
 *     provided with the distribution.                              *
 *  3. Neither the name of the copyright holder nor the names of    *
 *     its contributors may be used to endorse or promote products  *
 *     derived from this software without specific                  *
 *     prior written permission.                                    *
 *                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND           *
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,      *
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF         *
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE         *
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR            *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,     *
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,         *
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER *
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,      *
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)    *
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF      *
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
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
