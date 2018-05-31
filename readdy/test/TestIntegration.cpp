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
#include <readdy/common/numeric.h>

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
    EXPECT_LE(result.second/result.first, std::numeric_limits<readdy::scalar>::epsilon())
                        << "relative error should be <= machine precision ";
}

TEST(TestIntegration, IntegrateAdaptivePeriodicToMachinePrecision) {
    // integral in range [-a, +a] of odd function is zero, sin=odd, x^4=even, sin(x)*x^4=odd
    auto integrand = [](const readdy::scalar x) { return std::sin(x) * std::pow(x, 4); };
    auto result = readdy::util::integration::integrateAdaptive(integrand, -30., 30.,
                                                               std::numeric_limits<readdy::scalar>::epsilon(),
                                                               100);
    auto trueIntegral = 0.;
    EXPECT_FLOAT_EQ(result.first, trueIntegral);
    EXPECT_LE(result.second, std::numeric_limits<readdy::scalar>::epsilon())<< "error should be <= machine precision ";
}

TEST(TestIntegration, LimitsAreEqual) {
    auto integrand = [](const readdy::scalar x){return 10.;};
    auto result = readdy::util::integration::integrateAdaptive(integrand, 1., 1.);
    EXPECT_FLOAT_EQ(result.first, 0.);
    EXPECT_FLOAT_EQ(result.second, 0.);
}

TEST(TestIntegration, UnorderedLimits) {
    auto integrand = [](const readdy::scalar x){return 10.;};
    EXPECT_ANY_THROW(readdy::util::integration::integrateAdaptive(integrand, 1., 0.));
}

TEST(TestIntegration, RadialIntegralHarmonicRepulsion) {
    const readdy::scalar interactionDistance = 2.;
    const auto harmonicRepulsion = [&interactionDistance](const readdy::scalar x){
        if (x < 0.) {
            throw std::invalid_argument("only positive");
        }
        if (x < interactionDistance) {
            return 0.5 * 10. * std::pow(x - interactionDistance, 2.);
        } else {
            return 0.;
        }
    };
    auto integrand = [&harmonicRepulsion, &interactionDistance](const readdy::scalar x){
        return 4. * readdy::util::numeric::pi<readdy::scalar>() * x * x * std::exp(- harmonicRepulsion(x));
    };
    const readdy::scalar desiredRelativeError = 1e-12;
    auto result = readdy::util::integration::integrateAdaptive(integrand, 0., 5., desiredRelativeError);
    EXPECT_LE(result.second/result.first, desiredRelativeError);
}

TEST(TestIntegration, VaryingPeriodicity) {
    auto integrand = [](const readdy::scalar x) {return x * std::sin(1./x);};
    readdy::scalar desiredRelativeError = std::numeric_limits<readdy::scalar>::epsilon() * 1000;
    auto result = readdy::util::integration::integrateAdaptive(
            integrand, 0., readdy::util::numeric::pi<readdy::scalar>(), desiredRelativeError, 2000);
    readdy::scalar trueIntegral = 2.409156679836905734596964251490166243591343837595845857659;
    EXPECT_NEAR(result.first, trueIntegral, desiredRelativeError*trueIntegral);
    EXPECT_LE(result.second/result.first, desiredRelativeError);
}
