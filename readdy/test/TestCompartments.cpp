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
