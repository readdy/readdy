/********************************************************************
 * Copyright © 2018 Computational Molecular Biology Group,          *
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
 * « detailed description »
 *
 * @file TestIntegration.cpp
 * @brief « brief description »
 * @author chrisfroe
 * @date 24.05.18
 */

#include <gtest/gtest.h>
#include <readdy/common/integration.h>

TEST(TestIntegration, IntegratePolynomialExact) {
    // the used integration rules should yield exact results for integrating polynomials of order less than 2*201 - 1
    auto polynomial = [](const readdy::scalar x) { return  2.*x - 3.*std::pow(x, 2.) + 6.*std::pow(x, 5.); };
    const auto result = readdy::util::integration::integrate(polynomial, 0, 10);
    const auto numericIntegral = result.first;
    const auto errorEstimate = result.second;
    const readdy::scalar trueIntegral = 1e2 - 1e3 + 1e6;
    EXPECT_FLOAT_EQ(numericIntegral, trueIntegral) << "numeric result should be exact";
    EXPECT_FLOAT_EQ(errorEstimate, 0.) << "estimated error should be zero";
}

TEST(TestIntegration, IntegrateExponentialWithinEstimatedError) {
    auto integrand = [](const readdy::scalar x) { return std::exp(-x); };
    const auto result = readdy::util::integration::integrate(integrand, 0, 1);
    const auto numericIntegral = result.first;
    const auto errorEstimate = result.second;
    const auto trueIntegral = 1.0 - std::exp(-1.0);
    EXPECT_GT(errorEstimate, 0.) << "exponential cannot be integrated exactly";
    EXPECT_TRUE(std::abs(numericIntegral - trueIntegral) < errorEstimate);
}

TEST(TestIntegration, IntegrateAdaptiveExponentialToMachinePrecision) {
    throw nullptr;
}