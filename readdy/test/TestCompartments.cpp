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
 * @file TestCompartments.cpp
 * @brief Test implementation-independent execution of program 'Compartments'
 * @author chrisfroe
 * @date 18.10.16
 */

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <readdy/testing/Utils.h>
#include <readdy/testing/KernelTest.h>
#include <readdy/model/compartments/Compartments.h>

namespace m = readdy::model;

namespace {

class TestCompartments : public KernelTest {

};

TEST_P(TestCompartments, OneCompartmentOneConversionOneParticle) {
    auto &ctx = kernel->context();
    ctx.particleTypes().add("A", 1.);
    ctx.particleTypes().add("B", 1.);
    kernel->addParticle("A", readdy::Vec3(1, 0, 2));

    std::unordered_map<std::string, std::string> conversionsMap = {{"A", "B"}};
    ctx.compartments().addSphere(conversionsMap, "kugelrund", readdy::Vec3(0,0,0), 10., false);

    std::vector<std::string> typesToCount = {"A", "B"};
    auto &&obs = kernel->observe().nParticles(1, typesToCount);
    obs->evaluate();
    const auto &resultBefore = obs->getResult();
    EXPECT_THAT(resultBefore, ::testing::ElementsAre(1, 0)) << "Expect one A particle before program execution";

    auto &&evaluateCompartments = kernel->actions().evaluateCompartments();
    evaluateCompartments->perform();

    obs->evaluate();
    const auto &resultAfter = obs->getResult();
    EXPECT_THAT(resultAfter, ::testing::ElementsAre(0, 1)) << "Expect zero A particle after program execution";
}

TEST_P(TestCompartments, TwoCompartments) {
    // two compartments, four species A,B,C and D, in the end there should only be C and D particles
    auto &ctx = kernel->context();
    ctx.boxSize() = {{10, 10, 10}};
    ctx.particleTypes().add("A", 1.);
    ctx.particleTypes().add("B", 1.);
    ctx.particleTypes().add("C", 1.);
    ctx.particleTypes().add("D", 1.);
    auto &&comp = kernel->actions().evaluateCompartments();

    std::unordered_map<std::string, std::string> conversionsXPos = {{"A", "C"}, {"B", "C"}};
    std::unordered_map<std::string, std::string> conversionsXNeg = {{"A", "D"}, {"B", "D"}};
    ctx.compartments().addPlane(conversionsXPos, "XPos", readdy::Vec3(1,0,0), 0, true);
    ctx.compartments().addPlane(conversionsXNeg, "XNeg", readdy::Vec3(-1,0,0), 0, true);

    for (auto i = 0; i < 100; ++i) {
        kernel->addParticle("A", readdy::model::rnd::normal3<readdy::scalar>());
        kernel->addParticle("B", readdy::model::rnd::normal3<readdy::scalar>());
    }

    comp->perform();

    std::vector<std::string> typesToCount = {"A", "B", "C", "D"};
    auto &&obs = kernel->observe().nParticles(1, typesToCount);
    obs->evaluate();
    const auto &result = obs->getResult();
    EXPECT_EQ(result[0], 0) << "Expect no As";
    EXPECT_EQ(result[1], 0) << "Expect no Bs";

    const std::vector<std::string> typesC = {"C"};
    const std::vector<std::string> typesD = {"D"};
    auto &&posC = kernel->observe().positions(1, typesC);
    auto &&posD = kernel->observe().positions(1, typesD);
    posC->evaluate();
    posD->evaluate();
    const auto &resultC = posC->getResult();
    const auto &resultD = posD->getResult();
    for (auto pos : resultC) {
        EXPECT_GE(pos[0], 0) << "x position of Cs must be > 0";
    }
    for (auto pos : resultD) {
        EXPECT_LT(pos[0], 0) << "x position of Ds must be < 0";
    }
}

INSTANTIATE_TEST_CASE_P(TestCompartments, TestCompartments, ::testing::ValuesIn(readdy::testing::getKernelsToTest()));

}
