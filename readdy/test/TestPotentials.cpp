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
    kernel.context().particleTypes().add("A", static_cast<readdy::scalar>(1.));
    kernel.context().particleTypes().add("B", static_cast<readdy::scalar>(0.1));
    kernel.context().periodicBoundaryConditions() = {{false, false, false}};
    const unsigned int nParticlesA = 10;
    const unsigned int nParticlesB = 10;

    // add particles close to origin but not all exactly on 0,0,0
    for (auto i = 0; i < nParticlesA; ++i) kernel.addParticle("A", readdy::model::rnd::normal3<readdy::scalar>());
    for (auto i = 0; i < nParticlesB; ++i) kernel.addParticle("B", readdy::model::rnd::normal3<readdy::scalar>());
}

void run(readdy::model::Kernel &kernel, readdy::scalar timeStep) {
    unsigned int nSteps = 200;
    auto &&integrator = kernel.actions().eulerBDIntegrator(timeStep);
    auto &&nl = kernel.actions().updateNeighborList();
    auto &&forces = kernel.actions().calculateForces();
    nl->perform();
    for (readdy::time_step_type &&t = 0; t < nSteps; ++t) {
        forces->perform();
        integrator->perform();
        nl->perform();
    }
}

TEST_P(TestPotentials, TestParticlesStayInBox) {
    kernel->context().boxSize() = {{5, 5, 5}};
    const readdy::scalar timeStep = .005;

    setupParticles(*kernel);

    kernel->context().potentials().addHarmonicRepulsion("A", "B", .01, .1 + .01);

    std::array<std::string, 2> types{{"A", "B"}};
    //std::array<readdy::scalar, 2> radii {{.1, .01}};
    for (const auto &t : types) {
        // create cube potential that is spanned from (-1,-1,-1) to (1, 1, 1)
        readdy::Vec3 origin {-1, -1, -1};
        readdy::Vec3 extent {2, 2, 2};
        kernel->context().potentials().addBox(t, 10, origin, extent);
    }

    auto ppObs = kernel->observe().positions(1);
    readdy::Vec3 lowerBound{static_cast<readdy::scalar>(-2.5), static_cast<readdy::scalar>(-2.5),
                                   static_cast<readdy::scalar>(-2.5)},
            upperBound{2.5, 2.5, 2.5};
    ppObs->callback() = [lowerBound, upperBound](readdy::model::observables::Positions::result_type currentResult) {
        readdy::Vec3 avg{0, 0, 0};
        bool allWithinBounds = true;
        for (auto &&v : currentResult) {
            allWithinBounds &= v >= lowerBound && v <= upperBound;
            avg += v;
        }
        avg /= currentResult.size();
        if (!allWithinBounds) readdy::log::debug("Average position: {}", avg);
        ASSERT_TRUE(allWithinBounds);
    };
    auto connection = kernel->connectObservable(ppObs.get());
    run(*kernel, timeStep);
}

TEST_P(TestPotentials, TestParticleStayInSphere) {
    kernel->context().boxSize() = {{10, 10, 10}};
    const readdy::scalar timeStep = .005;

    setupParticles(*kernel);

    std::array<std::string, 2> types{{"A", "B"}};
    for (const auto &t : types) {
        readdy::Vec3 origin (0, 0, 0);
        kernel->context().potentials().addSphereIn(t, 20, origin, 3);
    }
    auto ppObs = kernel->observe().positions(1);
    const readdy::scalar maxDistFromOrigin = 4.0; // at kbt=1 and force_const=20 the RMSD in a well potential would be ~0.2
    const readdy::scalar maxDistFromOriginSquared = maxDistFromOrigin * maxDistFromOrigin;
    ppObs->callback() = [maxDistFromOriginSquared](readdy::model::observables::Positions::result_type currentResult) {
        readdy::Vec3 avg{0, 0, 0};
        bool allWithinBounds = true;
        for (auto &&v : currentResult) {
            const readdy::scalar distanceFromOriginSquared = v * v;
            allWithinBounds &= distanceFromOriginSquared < maxDistFromOriginSquared;
            avg += v;
        }
        avg /= currentResult.size();
        if (!allWithinBounds) readdy::log::debug("Average position: {}", avg);
        ASSERT_TRUE(allWithinBounds);
    };
    auto connection = kernel->connectObservable(ppObs.get());
    run(*kernel, timeStep);
}

TEST_P(TestPotentials, TestLennardJonesRepellent) {
    auto calculateForces = kernel->actions().calculateForces();
    // test system where the particles are closer together than they should be, i.e., the force should be repellent
    auto& ctx = kernel->context();
    // one particle type A
    ctx.particleTypes().add("A", 1.0);
    // large enough box
    ctx.boxSize() = {{10, 10, 10}};
    // particles are aligned in the x-y plane and have a distance of .09
    auto id0 = kernel->addParticle("A", {0, 0, 0});
    auto id1 = kernel->addParticle("A", {0, 0, .09});

    auto id2 = kernel->addParticle("A", {2, 0, 0});
    auto id3 = kernel->addParticle("A", {2, 0, .09});

    // the potential has exponents 3 and 2, a cutoff distance of 1.0, does not shift the energy, a well depth
    // of 1.0 and a zero-interaction distance of 0.1 (particle distance < sigma ==> repellent)
    kernel->context().potentials().addLennardJones("A", "A", 3, 2, 1.0, false, 1.0, .1);

    // record ids
    auto pObs = kernel->observe().particles(1);
    std::vector<readdy::model::Particle::id_type> ids;
    pObs->callback() = [&ids](const readdy::model::observables::Particles::result_type& result) {
        const auto& recordedIds = std::get<1>(result);
        ids.insert(ids.end(), recordedIds.begin(), recordedIds.end());
    };
    auto connParticles = kernel->connectObservable(pObs.get());
    // also record forces
    auto fObs = kernel->observe().forces(1);
    std::vector<readdy::Vec3> collectedForces;
    fObs->callback() = [&collectedForces](const readdy::model::observables::Forces::result_type& result) {
        collectedForces.insert(collectedForces.end(), result.begin(), result.end());
    };
    auto conn = kernel->connectObservable(fObs.get());

    // we need to update the neighbor list as this is a pair potential
    auto neighborList = kernel->actions().updateNeighborList(readdy::model::actions::UpdateNeighborList::init);
    neighborList->perform();

    auto updateNeighborList = kernel->actions().updateNeighborList(readdy::model::actions::UpdateNeighborList::update);
    updateNeighborList->perform();
    // calc forces
    calculateForces->perform();
    // give me results
    kernel->evaluateObservables(1);

    // the reference values were calculated numerically
    EXPECT_NEAR(kernel->stateModel().energy(), static_cast<readdy::scalar>(2.0 * 0.925925925926), 1e-6);
    auto id0Idx = std::find(ids.begin(), ids.end(), id0) - ids.begin();
    auto id1Idx = std::find(ids.begin(), ids.end(), id1) - ids.begin();
    auto id2Idx = std::find(ids.begin(), ids.end(), id2) - ids.begin();
    auto id3Idx = std::find(ids.begin(), ids.end(), id3) - ids.begin();
    readdy::Vec3 forceOnParticle0 {0, 0, static_cast<readdy::scalar>(-123.45679012)};
    readdy::Vec3 forceOnParticle1 {0, 0, static_cast<readdy::scalar>(123.45679012)};
    readdy::Vec3 forceOnParticle2 {0, 0, static_cast<readdy::scalar>(-123.45679012)};
    readdy::Vec3 forceOnParticle3 {0, 0, static_cast<readdy::scalar>(123.45679012)};

    if(kernel->singlePrecision()) {
        EXPECT_VEC3_NEAR(collectedForces[id0Idx], forceOnParticle0, 1e-4);
        EXPECT_VEC3_NEAR(collectedForces[id1Idx], forceOnParticle1, 1e-4);
        EXPECT_VEC3_NEAR(collectedForces[id2Idx], forceOnParticle2, 1e-4);
        EXPECT_VEC3_NEAR(collectedForces[id3Idx], forceOnParticle3, 1e-4);
    } else {
        EXPECT_VEC3_NEAR(collectedForces[id0Idx], forceOnParticle0, 1e-6);
        EXPECT_VEC3_NEAR(collectedForces[id1Idx], forceOnParticle1, 1e-6);
        EXPECT_VEC3_NEAR(collectedForces[id2Idx], forceOnParticle2, 1e-6);
        EXPECT_VEC3_NEAR(collectedForces[id3Idx], forceOnParticle3, 1e-6);
    }
}

TEST_P(TestPotentials, ScreenedElectrostatics) {
    auto calculateForces = kernel->actions().calculateForces();
    auto &ctx = kernel->context();
    ctx.periodicBoundaryConditions() = {{false, false, false}};
    ctx.particleTypes().add("A", 1.0);
    ctx.boxSize() = {{10, 10, 10}};
    ctx.potentials().addBox("A", .001, {-4.9, -4.9, -4.9}, {9.8, 9.8, 9.8});
    // distance of particles is 2.56515106768
    auto id0 = kernel->addParticle("A", {0, 0, 0});
    auto id1 = kernel->addParticle("A", {1.2, 1.5, -1.7});
    auto electrostaticStrength = static_cast<readdy::scalar>(-1.);
    auto screeningDepth = static_cast<readdy::scalar>(1.);
    auto repulsionStrength = static_cast<readdy::scalar>(1.);
    auto sigma = static_cast<readdy::scalar>(1.);
    auto cutoff = static_cast<readdy::scalar>(8.);
    unsigned int exponent = 6;
    kernel->context().potentials().addScreenedElectrostatics("A", "A", electrostaticStrength, 1. / screeningDepth,
                                                                                 repulsionStrength, sigma, exponent, cutoff);
    // record ids to get data-structure-indexes of the two particles later on
    auto pObs = kernel->observe().particles(1);
    std::vector<readdy::model::Particle::id_type> ids;
    pObs->callback() = [&ids](const readdy::model::observables::Particles::result_type &result) {
        const auto &recordedIds = std::get<1>(result);
        ids.insert(ids.end(), recordedIds.begin(), recordedIds.end());
    };
    auto connParticles = kernel->connectObservable(pObs.get());
    // also record forces
    auto fObs = kernel->observe().forces(1);
    std::vector<readdy::Vec3> collectedForces;
    fObs->callback() = [&collectedForces](const readdy::model::observables::Forces::result_type &result) {
        collectedForces.insert(collectedForces.end(), result.begin(), result.end());
    };
    auto conn = kernel->connectObservable(fObs.get());

    // we need to update the neighbor list as this is a pair potential
    auto neighborList = kernel->actions().updateNeighborList();
    neighborList->perform();
    // calc forces
    calculateForces->perform();
    // give me results
    kernel->evaluateObservables(1);

    // the reference values were calculated numerically
    if(kernel->singlePrecision()) {
        EXPECT_FLOAT_EQ(kernel->stateModel().energy(), static_cast<readdy::scalar>(-0.0264715664281));
    } else {
        EXPECT_FLOAT_EQ(kernel->stateModel().energy(), static_cast<readdy::scalar>(-0.0264715664281));
    }
    ptrdiff_t id0Idx = std::find(ids.begin(), ids.end(), id0) - ids.begin();
    ptrdiff_t id1Idx = std::find(ids.begin(), ids.end(), id1) - ids.begin();
    readdy::Vec3 forceOnParticle0{0.01565262, 0.01956577, static_cast<readdy::scalar>(-0.02217454)};
    readdy::Vec3 forceOnParticle1{static_cast<readdy::scalar>(-0.01565262),
                                         static_cast<readdy::scalar>(-0.01956577), 0.02217454};
    EXPECT_VEC3_NEAR(collectedForces[id0Idx], forceOnParticle0, 1e-8);
    EXPECT_VEC3_NEAR(collectedForces[id1Idx], forceOnParticle1, 1e-8);
}

TEST_P(TestPotentials, SphericalMembrane) {
    // Combine SphereIn and SphereOut to build a 2D spherical manifold
    auto &ctx = kernel->context();
    ctx.particleTypes().add("A", 1.0);
    ctx.boxSize() =  {{10, 10, 10}};
    // add two particles, one outside, one inside the sphere
    auto id0 = kernel->addParticle("A", {2., 1., 1.});
    auto id1 = kernel->addParticle("A", {4., 3., -3.});
    auto forceConstant = static_cast<readdy::scalar>(1.);
    auto radius = static_cast<readdy::scalar>(3.);
    readdy::Vec3 origin = {static_cast<readdy::scalar>(1.), static_cast<readdy::scalar>(0.),
                                  static_cast<readdy::scalar>(0.)};
    kernel->context().potentials().addSphereOut("A", forceConstant, origin, radius);
    kernel->context().potentials().addSphereIn("A", forceConstant, origin, radius);
    // record ids to get data-structure-indexes of the two particles later on
    auto pObs = kernel->observe().particles(1);
    std::vector<readdy::model::Particle::id_type> ids;
    pObs->callback() = [&ids](const readdy::model::observables::Particles::result_type &result) {
        const auto &recordedIds = std::get<1>(result);
        ids.insert(ids.end(), recordedIds.begin(), recordedIds.end());
    };
    auto connParticles = kernel->connectObservable(pObs.get());
    // also record forces
    auto fObs = kernel->observe().forces(1);
    std::vector<readdy::Vec3> collectedForces;
    fObs->callback() = [&collectedForces](const readdy::model::observables::Forces::result_type &result) {
        collectedForces.insert(collectedForces.end(), result.begin(), result.end());
    };
    auto conn = kernel->connectObservable(fObs.get());

    // we need to update the neighbor list as this is a pair potential
    auto neighborList = kernel->actions().updateNeighborList();
    neighborList->perform();
    auto updateNeighborList = kernel->actions().updateNeighborList(readdy::model::actions::UpdateNeighborList::Operation::update);
    updateNeighborList->perform();
    // calc forces
    auto calculateForces = kernel->actions().calculateForces();
    calculateForces->perform();
    // give me results
    kernel->evaluateObservables(1);

    // the reference values were calculated numerically
    ASSERT_FLOAT_EQ(kernel->stateModel().energy(), 0.803847577293 + 2.41154273188);
    ptrdiff_t id0Idx = std::find(ids.begin(), ids.end(), id0) - ids.begin();
    ptrdiff_t id1Idx = std::find(ids.begin(), ids.end(), id1) - ids.begin();
    readdy::Vec3 forceOnParticle0{static_cast<readdy::scalar>(0.73205081),
                                         static_cast<readdy::scalar>(0.73205081),
                                         static_cast<readdy::scalar>(0.73205081)};
    readdy::Vec3 forceOnParticle1{static_cast<readdy::scalar>(-1.26794919),
                                         static_cast<readdy::scalar>(-1.26794919),
                                         static_cast<readdy::scalar>(1.26794919)};
    EXPECT_VEC3_NEAR(collectedForces[id0Idx], forceOnParticle0, 1e-6);
    EXPECT_VEC3_NEAR(collectedForces[id1Idx], forceOnParticle1, 1e-6);
}

TEST_P(TestPotentials, SphericalBarrier) {
    auto &ctx = kernel->context();
    auto calculateForces = kernel->actions().calculateForces();
    ctx.particleTypes().add("A", 1.0);
    ctx.boxSize() = {{10, 10, 10}};
    // add two particles, one on the outer edge getting pushed outside, one inside the sphere unaffected
    auto id0 = kernel->addParticle("A", {2.1, 1., 1.});
    auto id1 = kernel->addParticle("A", {1.1, 1., 1.});
    readdy::Vec3 origin = {1.0, 0.9, 1.0};
    double radius = 1.;
    double height = 2.;
    double width = 0.3;
    kernel->context().potentials().addSphericalBarrier("A", height, width, origin, radius);
    // record ids to get data-structure-indexes of the two particles later on
    auto pObs = kernel->observe().particles(1);
    std::vector<readdy::model::Particle::id_type> ids;
    pObs->callback() = [&ids](const readdy::model::observables::Particles::result_type &result) {
        const auto &recordedIds = std::get<1>(result);
        ids.insert(ids.end(), recordedIds.begin(), recordedIds.end());
    };
    auto connParticles = kernel->connectObservable(pObs.get());
    // also record forces
    auto fObs = kernel->observe().forces(1);
    std::vector<readdy::Vec3> collectedForces;
    fObs->callback() = [&collectedForces](const readdy::model::observables::Forces::result_type &result) {
        collectedForces.insert(collectedForces.end(), result.begin(), result.end());
    };
    auto connForces = kernel->connectObservable(fObs.get());

    kernel->initialize();

    calculateForces->perform();
    kernel->evaluateObservables(1);

    // the reference values were calculated numerically
    ASSERT_FLOAT_EQ(kernel->stateModel().energy(), 1.51432015278);
    auto id0Idx = std::find(ids.begin(), ids.end(), id0) - ids.begin();
    auto id1Idx = std::find(ids.begin(), ids.end(), id1) - ids.begin();
    readdy::Vec3 forceOnParticle0{static_cast<readdy::scalar>(9.2539372), static_cast<readdy::scalar>(0.84126702),
                                         static_cast<readdy::scalar>(0.)};
    readdy::Vec3 forceOnParticle1{static_cast<readdy::scalar>(0.), static_cast<readdy::scalar>(0.),
                                         static_cast<readdy::scalar>(0.)};
    EXPECT_VEC3_NEAR(collectedForces[id0Idx], forceOnParticle0, kernel->doublePrecision() ? 1e-8 : 1e-5);
    EXPECT_VEC3_NEAR(collectedForces[id1Idx], forceOnParticle1, kernel->doublePrecision() ? 1e-8 : 1e-5);
}

INSTANTIATE_TEST_CASE_P(TestPotentials, TestPotentials,
                        ::testing::ValuesIn(readdy::testing::getKernelsToTest()));
}
