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
    for (readdy::time_step_type &&t = 0; t < nSteps; ++t) {
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
        kernel->registerPotential<readdy::model::potentials::Cube>(t, 10, origin, extent, true);
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
        if (!allWithinBounds) readdy::log::debug("Average position: {}", avg);
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
        kernel->registerPotential<readdy::model::potentials::SphereIn>(t, 20, origin, 3);
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
        if (!allWithinBounds) readdy::log::debug("Average position: {}", avg);
        ASSERT_TRUE(allWithinBounds);
    });
    auto connection = kernel->connectObservable(ppObs.get());
    run(*kernel, timeStep);
}

TEST_P(TestPotentials, TestLennardJonesRepellent) {
    // test system where the particles are closer together than they should be, i.e., the force should be repellent
    auto& ctx = kernel->getKernelContext();
    // one particle type A
    ctx.registerParticleType("A", 1.0, 1.0);
    // large enough box
    ctx.setBoxSize(10, 10, 10);
    // particles are aligned in the x-y plane and have a distance of .09
    auto id0 = kernel->addParticle("A", {0, 0, 0});
    auto id1 = kernel->addParticle("A", {0, 0, .09});

    // the potential has exponents 3 and 2, a cutoff distance of 1.0, does not shift the energy, a well depth
    // of 1.0 and a zero-interaction distance of 0.1 (particle distance < sigma ==> repellent)
    kernel->registerPotential<readdy::model::potentials::LennardJones>("A", "A", 3, 2, 1.0, false, 1.0, .1);

    // record ids
    auto pObs = kernel->createObservable<readdy::model::observables::Particles>(1);
    std::vector<readdy::model::Particle::id_type> ids;
    pObs->setCallback([&ids](const readdy::model::observables::Particles::result_t& result) {
        const auto& recordedIds = std::get<1>(result);
        ids.insert(ids.end(), recordedIds.begin(), recordedIds.end());
    });
    auto connParticles = kernel->connectObservable(pObs.get());
    // also record forces
    auto fObs = kernel->createObservable<readdy::model::observables::Forces>(1);
    std::vector<readdy::model::Vec3> collectedForces;
    fObs->setCallback([&collectedForces](const readdy::model::observables::Forces::result_t& result) {
        collectedForces.insert(collectedForces.end(), result.begin(), result.end());
    });
    auto conn = kernel->connectObservable(fObs.get());
    // configure the context
    ctx.configure();

    // we need to update the neighbor list as this is a pair potential
    auto neighborList = kernel->createAction<readdy::model::actions::UpdateNeighborList>();
    neighborList->perform();
    // calc forces
    kernel->getKernelStateModel().calculateForces();
    // give me results
    kernel->evaluateObservables(1);

    // the reference values were calculated numerically
    ASSERT_FLOAT_EQ(kernel->getKernelStateModel().getEnergy(), 0.925925925926);
    ptrdiff_t id0Idx = std::find(ids.begin(), ids.end(), id0) - ids.begin();
    ptrdiff_t id1Idx = std::find(ids.begin(), ids.end(), id1) - ids.begin();
    readdy::model::Vec3 forceOnParticle0 {0, 0, -123.45679012};
    readdy::model::Vec3 forceOnParticle1 {0, 0, 123.45679012};
    EXPECT_VEC3_NEAR(collectedForces[id0Idx], forceOnParticle0, 1e-6);
    EXPECT_VEC3_NEAR(collectedForces[id1Idx], forceOnParticle1, 1e-6);
}

TEST_P(TestPotentials, ShieldedElectrostatics) {
    auto &ctx = kernel->getKernelContext();
    ctx.registerParticleType("A", 1.0, 1.0);
    ctx.setBoxSize(10, 10, 10);
    // distance of particles is 2.56515106768
    auto id0 = kernel->addParticle("A", {0, 0, 0});
    auto id1 = kernel->addParticle("A", {1.2, 1.5, -1.7});
    double electrostaticStrength = -1.;
    double screeningDepth = 1.;
    double repulsionStrength = 1.;
    double sigma = 1.;
    double cutoff = 8.;
    unsigned int exponent = 6;
    kernel->registerPotential<readdy::model::potentials::ShieldedElectrostatics>("A", "A", electrostaticStrength, 1. / screeningDepth,
                                                                                 repulsionStrength, sigma, exponent, cutoff);
    // record ids to get data-structure-indexes of the two particles later on
    auto pObs = kernel->createObservable<readdy::model::observables::Particles>(1);
    std::vector<readdy::model::Particle::id_type> ids;
    pObs->setCallback([&ids](const readdy::model::observables::Particles::result_t &result) {
        const auto &recordedIds = std::get<1>(result);
        ids.insert(ids.end(), recordedIds.begin(), recordedIds.end());
    });
    auto connParticles = kernel->connectObservable(pObs.get());
    // also record forces
    auto fObs = kernel->createObservable<readdy::model::observables::Forces>(1);
    std::vector<readdy::model::Vec3> collectedForces;
    fObs->setCallback([&collectedForces](const readdy::model::observables::Forces::result_t &result) {
        collectedForces.insert(collectedForces.end(), result.begin(), result.end());
    });
    auto conn = kernel->connectObservable(fObs.get());
    // configure the context
    ctx.configure();

    // we need to update the neighbor list as this is a pair potential
    auto neighborList = kernel->createAction<readdy::model::actions::UpdateNeighborList>();
    neighborList->perform();
    // calc forces
    kernel->getKernelStateModel().calculateForces();
    // give me results
    kernel->evaluateObservables(1);

    // the reference values were calculated numerically
    ASSERT_FLOAT_EQ(kernel->getKernelStateModel().getEnergy(), -0.0264715664281);
    ptrdiff_t id0Idx = std::find(ids.begin(), ids.end(), id0) - ids.begin();
    ptrdiff_t id1Idx = std::find(ids.begin(), ids.end(), id1) - ids.begin();
    readdy::model::Vec3 forceOnParticle0{0.01565262, 0.01956577, -0.02217454};
    readdy::model::Vec3 forceOnParticle1{-0.01565262, -0.01956577, 0.02217454};
    EXPECT_VEC3_NEAR(collectedForces[id0Idx], forceOnParticle0, 1e-8);
    EXPECT_VEC3_NEAR(collectedForces[id1Idx], forceOnParticle1, 1e-8);
}

TEST_P(TestPotentials, SphericalMembrane) {
    // Combine SphereIn and SphereOut to build a 2D spherical manifold
    auto &ctx = kernel->getKernelContext();
    ctx.registerParticleType("A", 1.0, 1.0);
    ctx.setBoxSize(10, 10, 10);
    // add two particles, one outside, one inside the sphere
    auto id0 = kernel->addParticle("A", {2., 1., 1.});
    auto id1 = kernel->addParticle("A", {4., 3., -3.});
    double forceConstant = 1.;
    double radius = 3.;
    readdy::model::Vec3 origin = {1.,0.,0.};
    kernel->registerPotential<readdy::model::potentials::SphereOut>("A", forceConstant, origin, radius);
    kernel->registerPotential<readdy::model::potentials::SphereIn>("A", forceConstant, origin, radius);
    // record ids to get data-structure-indexes of the two particles later on
    auto pObs = kernel->createObservable<readdy::model::observables::Particles>(1);
    std::vector<readdy::model::Particle::id_type> ids;
    pObs->setCallback([&ids](const readdy::model::observables::Particles::result_t &result) {
        const auto &recordedIds = std::get<1>(result);
        ids.insert(ids.end(), recordedIds.begin(), recordedIds.end());
    });
    auto connParticles = kernel->connectObservable(pObs.get());
    // also record forces
    auto fObs = kernel->createObservable<readdy::model::observables::Forces>(1);
    std::vector<readdy::model::Vec3> collectedForces;
    fObs->setCallback([&collectedForces](const readdy::model::observables::Forces::result_t &result) {
        collectedForces.insert(collectedForces.end(), result.begin(), result.end());
    });
    auto conn = kernel->connectObservable(fObs.get());
    // configure the context
    ctx.configure();

    // we need to update the neighbor list as this is a pair potential
    auto neighborList = kernel->createAction<readdy::model::actions::UpdateNeighborList>();
    neighborList->perform();
    // calc forces
    kernel->getKernelStateModel().calculateForces();
    // give me results
    kernel->evaluateObservables(1);

    // the reference values were calculated numerically
    ASSERT_FLOAT_EQ(kernel->getKernelStateModel().getEnergy(), 0.803847577293 + 2.41154273188);
    ptrdiff_t id0Idx = std::find(ids.begin(), ids.end(), id0) - ids.begin();
    ptrdiff_t id1Idx = std::find(ids.begin(), ids.end(), id1) - ids.begin();
    readdy::model::Vec3 forceOnParticle0{0.73205081, 0.73205081, 0.73205081};
    readdy::model::Vec3 forceOnParticle1{-1.26794919, -1.26794919, 1.26794919};
    EXPECT_VEC3_NEAR(collectedForces[id0Idx], forceOnParticle0, 1e-8);
    EXPECT_VEC3_NEAR(collectedForces[id1Idx], forceOnParticle1, 1e-8);
}

INSTANTIATE_TEST_CASE_P(TestPotentials, TestPotentials,
                        ::testing::ValuesIn(readdy::testing::getKernelsToTest()));
}
