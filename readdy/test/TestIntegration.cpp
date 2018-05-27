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
 * @file TestIntegration.cpp
 * @brief Tests for numerical integration
 * @author chrisfroe
 * @date 24.05.18
 */

#include <gtest/gtest.h>
#include <readdy/common/common.h>
#include <readdy/common/integration.h>

TEST(TestIntegration, IntegratePolynomialExact) {
    // the used integration rules should yield exact results for integrating polynomials of order less than 2*201 - 1
    auto polynomial = [](const readdy::scalar x) { return 2. * x - 3. * std::pow(x, 2.) + 6. * std::pow(x, 5.); };
    auto result = readdy::util::integration::integrate(polynomial, 0., 10.);
    auto numericIntegral = result.first;
    auto errorEstimate = result.second;
    readdy::scalar trueIntegral = 1e2 - 1e3 + 1e6;
    EXPECT_FLOAT_EQ(numericIntegral, trueIntegral) << "numeric result should be exact";
    EXPECT_FLOAT_EQ(errorEstimate, 0.) << "estimated error should be zero";
}

TEST(TestIntegration, IntegrateExponentialWithinEstimatedError) {
    auto integrand = [](const readdy::scalar x) { return std::exp(-x); };
    auto result = readdy::util::integration::integrate(integrand, 0., 1.);
    auto numericIntegral = result.first;
    auto errorEstimate = result.second;
    auto trueIntegral = 1.0 - std::exp(-1.0);
    EXPECT_GT(errorEstimate, 0.) << "exponential cannot be integrated exactly";
    EXPECT_TRUE(std::abs(numericIntegral - trueIntegral) < errorEstimate);
}

TEST(TestIntegration, IntegrateAdaptiveExponentialToMachinePrecision) {
    auto integrand = [](const readdy::scalar x) { return std::exp(-x); };
    auto result = readdy::util::integration::integrateAdaptive(integrand, 0., 1.,
                                                               std::numeric_limits<readdy::scalar>::epsilon(),
                                                               100);
    auto trueIntegral = 1.0 - std::exp(-1.0);
    EXPECT_FLOAT_EQ(result.first, trueIntegral);
    EXPECT_LE(result.second, std::numeric_limits<readdy::scalar>::epsilon()) << "error should be <= machine precision ";
}

TEST(TestIntegration, IntegrateAdaptivePeriodicToMachinePrecision) {
    // integral in range [-a, +a] of odd function is zero, sin=odd, x^4=even, sin(x)*x^4=odd
    auto integrand = [](const readdy::scalar x) { return std::sin(x) * std::pow(x, 4); };
    auto result = readdy::util::integration::integrateAdaptive(integrand, -30., 30.,
                                                               std::numeric_limits<readdy::scalar>::epsilon(),
                                                               100);
    auto trueIntegral = 0.;
    EXPECT_FLOAT_EQ(result.first, trueIntegral);
    EXPECT_LE(result.second, std::numeric_limits<readdy::scalar>::epsilon()) << "error should be <= machine precision ";
}