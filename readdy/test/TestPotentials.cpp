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
 * Test potentials of first and second order via small particle systems. First order potentials
 * often assure that particles stay in a certain volume of the simulation box, which is checked for.
 *
 * @file TestPotentials.cpp
 * @brief Check correct behavior of particles subject to potentials.
 * @author clonker
 * @author chrisfroe
 * @date 27.06.16
 */

#include <gtest/gtest.h>
#include <readdy/plugin/KernelProvider.h>
#include <readdy/testing/KernelTest.h>
#include <readdy/testing/Utils.h>

namespace {

class TestPotentials : public KernelTest {

};

void setupParticles(readdy::model::Kernel &kernel) {
    kernel.getKernelContext().registerParticleType("A", 1., 0.1);
    kernel.getKernelContext().registerParticleType("B", 0.1, 0.01);
    kernel.getKernelContext().setPeriodicBoundary(false, false, false);
    const unsigned int nParticlesA = 10;
    const unsigned int nParticlesB = 10;

    // add particles close to origin but not all exactly on 0,0,0
    for (auto i = 0; i < nParticlesA; ++i) kernel.addParticle("A", readdy::model::rnd::normal3());
    for (auto i = 0; i < nParticlesB; ++i) kernel.addParticle("B", readdy::model::rnd::normal3());
}

void run(readdy::model::Kernel &kernel, double timeStep) {
    unsigned int nSteps = 2000;
    auto &&integrator = kernel.createAction<readdy::model::actions::EulerBDIntegrator>(timeStep);
    auto &&nl = kernel.createAction<readdy::model::actions::UpdateNeighborList>();
    auto &&forces = kernel.createAction<readdy::model::actions::CalculateForces>();
    nl->perform();
    for (readdy::model::observables::time_step_type &&t = 0; t < nSteps; ++t) {
        forces->perform();
        integrator->perform();
        nl->perform();
    }
}

TEST_P(TestPotentials, TestParticlesStayInBox) {
    kernel->getKernelContext().setBoxSize(5, 5, 5);
    const double timeStep = .005;

    setupParticles(*kernel);

    kernel->registerPotential<readdy::model::potentials::HarmonicRepulsion>("A", "B", .01);

    std::array<std::string, 2> types{{"A", "B"}};
    for (auto t : types) {
        // create cube potential that is spanned from (-1,-1,-1) to (1, 1, 1)
        readdy::model::Vec3 origin {-1, -1, -1};
        readdy::model::Vec3 extent {2, 2, 2};
        kernel->registerPotential<readdy::model::potentials::CubePotential>(t, 10, origin, extent, true);
    }

    auto ppObs = kernel->createObservable<readdy::model::observables::Positions>(1);
    readdy::model::Vec3 lowerBound{-2.5, -2.5, -2.5}, upperBound{2.5, 2.5, 2.5};
    ppObs->setCallback([lowerBound, upperBound](readdy::model::observables::Positions::result_t currentResult) {
        readdy::model::Vec3 avg{0, 0, 0};
        bool allWithinBounds = true;
        for (auto &&v : currentResult) {
            allWithinBounds &= v >= lowerBound && v <= upperBound;
            avg += v;
        }
        avg /= currentResult.size();
        if (!allWithinBounds) readdy::log::console()->debug("Average position: {}", avg);
        ASSERT_TRUE(allWithinBounds);
    });
    auto connection = kernel->connectObservable(ppObs.get());
    run(*kernel, timeStep);
}

TEST_P(TestPotentials, TestParticleStayInSphere) {
    kernel->getKernelContext().setBoxSize(10, 10, 10);
    const double timeStep = .005;

    setupParticles(*kernel);

    std::array<std::string, 2> types{{"A", "B"}};
    for (auto t : types) {
        readdy::model::Vec3 origin (0, 0, 0);
        kernel->registerPotential<readdy::model::potentials::SpherePotential>(t, 20, origin, 3);
    }
    auto ppObs = kernel->createObservable<readdy::model::observables::Positions>(1);
    const double maxDistFromOrigin = 4.0; // at kbt=1 and force_const=20 the RMSD in a well potential would be ~0.2
    const double maxDistFromOriginSquared = maxDistFromOrigin * maxDistFromOrigin;
    ppObs->setCallback([maxDistFromOriginSquared](readdy::model::observables::Positions::result_t currentResult) {
        readdy::model::Vec3 avg{0, 0, 0};
        bool allWithinBounds = true;
        for (auto &&v : currentResult) {
            const double distanceFromOriginSquared = v * v;
            allWithinBounds &= distanceFromOriginSquared < maxDistFromOriginSquared;
            avg += v;
        }
        avg /= currentResult.size();
        if (!allWithinBounds) readdy::log::console()->debug("Average position: {}", avg);
        ASSERT_TRUE(allWithinBounds);
    });
    auto connection = kernel->connectObservable(ppObs.get());
    run(*kernel, timeStep);
}

INSTANTIATE_TEST_CASE_P(TestPotentials, TestPotentials,
                        ::testing::ValuesIn(readdy::testing::getKernelsToTest()));
}
