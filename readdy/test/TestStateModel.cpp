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
 * @file TestStateModel.cpp
 * @brief Test the methods that manipulate time-dependent simulation data for kernels.
 * @author clonker
 * @author chrisfroe
 * @date 25.08.16
 * @todo check force calculation through periodic boundary
 */
#include <gtest/gtest.h>
#include <readdy/testing/Utils.h>
#include <readdy/testing/KernelTest.h>

namespace m = readdy::model;

namespace {

class TestStateModel : public KernelTest {

};

TEST_P(TestStateModel, CalculateForcesTwoParticles) {
    m::Context &ctx = kernel->context();
    auto &stateModel = kernel->stateModel();

    auto obs = kernel->observe().forces(1);
    auto conn = kernel->connectObservable(obs.get());
    // two A particles with radius 1. -> cutoff 2, distance 1.8 -> r-r_0 = 0.2 -> force = 0.2
    ctx.particleTypes().add("A", 1.0);
    ctx.boxSize() = {{4., 4., 4.}};
    ctx.periodicBoundaryConditions() = {{false, false, false}};

    ctx.potentials().addBox("A", .001, {-1.9, -1.9, -1.9}, {3.8, 3.8, 3.8});

    kernel->context().potentials().addHarmonicRepulsion("A", "A", 1.0, 2.0);
    kernel->context().potentials().addBox("A", .01, {-1.9, -1.9, -1.9}, {3.8, 3.8, 3.8});

    kernel->initialize();
    auto typeIdA = ctx.particleTypes().idOf("A");
    auto twoParticles = std::vector<m::Particle> {m::Particle(0., 0., 0., typeIdA), m::Particle(0., 0., 1.8, typeIdA)};

    stateModel.addParticles(twoParticles);
    stateModel.initializeNeighborList(0.);
    stateModel.updateNeighborList();

    auto calculateForces = kernel->actions().calculateForces();
    calculateForces->perform();
    calculateForces->perform(); // calculating twice should yield the same result. force and energy must not accumulate
    // check results
    obs->evaluate();
    auto forcesIt = obs->getResult().begin();
    if(kernel->doublePrecision()) {
        EXPECT_VEC3_EQ(*forcesIt, readdy::Vec3(0, 0, -0.2));
    } else {
        EXPECT_FVEC3_EQ(*forcesIt, readdy::Vec3(0, 0, -0.2));
    }
    ++forcesIt;
    if(kernel->doublePrecision()) {
        EXPECT_VEC3_EQ(*forcesIt, readdy::Vec3(0, 0, 0.2));
    } else {
        EXPECT_FVEC3_EQ(*forcesIt, readdy::Vec3(0, 0, 0.2));
    }
    if(kernel->doublePrecision()) {
        EXPECT_DOUBLE_EQ(stateModel.energy(), 0.02);
    } else {
        EXPECT_NEAR(stateModel.energy(), 0.02, 1e-8);
    }
}

TEST_P(TestStateModel, CalculateForcesRepulsion) {
    m::Context &ctx = kernel->context();
    auto calculateForces = kernel->actions().calculateForces();
    auto &stateModel = kernel->stateModel();

    // similar situation as before but now with repulsion between A and B
    auto obs = kernel->observe().forces(1);
    auto conn = kernel->connectObservable(obs.get());
    ctx.particleTypes().add("A", 1.0);
    ctx.particleTypes().add("B", 1.0);
    ctx.boxSize() = {{10., 10., 10.}};
    ctx.periodicBoundaryConditions() = {{true, true, false}};
    ctx.potentials().addHarmonicRepulsion("A", "B", 1.0, 3.0);
    kernel->context().potentials().addBox("A", .01, {-4.9, -4.9, -4.9}, {9.8, 9.8, 9.8});
    kernel->context().potentials().addBox("B", .01, {-4.9, -4.9, -4.9}, {9.8, 9.8, 9.8});

    auto typeIdA = ctx.particleTypes().idOf("A");
    auto typeIdB = ctx.particleTypes().idOf("B");
    /**
     * There are 6 particles. 0-2 are A particles. 3-5 are B particles.
     * The second B particle is a bit further away
     * -> There are forces between all AB pairs except with the second B particle, namely particle 4
     */
    auto particlesA = std::vector<m::Particle> {
            m::Particle(0, 0, 0, typeIdA), m::Particle(0, 0.8, 0, typeIdA), m::Particle(0.2, 0, -0.2, typeIdA)
    };
    auto particlesB = std::vector<m::Particle> {
            m::Particle(0, 0, 0.1, typeIdB), m::Particle(-1.8, -4.0, 1.8, typeIdB), m::Particle(0.5, 0, 0, typeIdB)
    };
    std::vector<m::Particle::id_type> ids;
    {
        std::for_each(particlesA.begin(), particlesA.end(), [&ids](const m::Particle& p) { ids.push_back(p.getId());});
        std::for_each(particlesB.begin(), particlesB.end(), [&ids](const m::Particle& p) { ids.push_back(p.getId());});
    }
    stateModel.addParticles(particlesA);
    stateModel.addParticles(particlesB);
    kernel->initialize();
    stateModel.initializeNeighborList(0);
    stateModel.updateNeighborList();
    calculateForces->perform();

    // by hand calculated expectations
    const readdy::scalar energy03 = 4.205;
    const readdy::scalar energy05 = 3.125;
    const readdy::scalar energy13 = 2.4063226755104354;
    const readdy::scalar energy15 = 2.1148056603830194;
    const readdy::scalar energy23 = 3.4833346173608031;
    const readdy::scalar energy25 = 3.4833346173608031;
    const readdy::Vec3 force03(0, 0, -2.9);
    const readdy::Vec3 force05(-2.5, 0, 0);
    const readdy::Vec3 force13(0, 2.1768336301410027, -0.27210420376762534);
    const readdy::Vec3 force15(-1.08999682000954, 1.743994912015264, 0);
    const readdy::Vec3 force23(1.4641005886756873, 0, -2.1961508830135306);
    const readdy::Vec3 force25(-2.1961508830135306, 0, -1.4641005886756873);

    // check results
    obs->evaluate();
    const auto particles = stateModel.getParticles();
    const auto& forces = obs->getResult();
    std::size_t idx = 0;
    for(const auto& particle : particles) {
        if(particle.getId() == ids.at(0)) {
            if(kernel->doublePrecision()) {
                EXPECT_VEC3_EQ(forces.at(idx), force03 + force05) << "force on particle 0 = force03 + force05";
            } else {
                EXPECT_FVEC3_EQ(forces.at(idx), force03 + force05) << "force on particle 0 = force03 + force05";
            }
        } else if(particle.getId() == ids.at(1)) {
            if(kernel->doublePrecision()) {
                EXPECT_VEC3_EQ(forces.at(idx), force13 + force15) << "force on particle 1 = force13 + force15";
            } else {
                EXPECT_FVEC3_EQ(forces.at(idx), force13 + force15) << "force on particle 1 = force13 + force15";
            }
        } else if(particle.getId() == ids.at(2)) {
            if(kernel->doublePrecision()) {
                EXPECT_VEC3_EQ(forces.at(idx), force23 + force25) << "force on particle 2 = force23 + force25";
            } else {
                EXPECT_FVEC3_EQ(forces.at(idx), force23 + force25) << "force on particle 2 = force23 + force25";
            }
        } else if(particle.getId() == ids.at(3)) {
            if(kernel->doublePrecision()) {
                EXPECT_VEC3_EQ(forces.at(idx), (-1. * force03) - force13 - force23)
                                    << "force on particle 3 = - force03 - force13 - force23";
            } else {
                EXPECT_FVEC3_EQ(forces.at(idx), (-1. * force03) - force13 - force23)
                                    << "force on particle 3 = - force03 - force13 - force23";
            }
        } else if(particle.getId() == ids.at(4)) {
            if(kernel->doublePrecision()) {
                EXPECT_VEC3_EQ(forces.at(idx), readdy::Vec3(0, 0, 0)) << "force on particle 4 = 0";
            } else {
                EXPECT_FVEC3_EQ(forces.at(idx), readdy::Vec3(0, 0, 0)) << "force on particle 4 = 0";
            }
        } else if(particle.getId() == ids.at(5)) {
            if(kernel->doublePrecision()) {
                EXPECT_VEC3_EQ(forces.at(idx), (-1. * force05) - force15 - force25)
                                    << "force on particle 5 = - force05 - force15 - force25";
            } else {
                EXPECT_FVEC3_EQ(forces.at(idx), (-1. * force05) - force15 - force25)
                                    << "force on particle 5 = - force05 - force15 - force25";
            }
        } else {
            readdy::log::error("Got an unexpected particle id: {}", particle.getId());
            FAIL() << "Got an unexpected particle id: " << particle.getId();
        }
        ++idx;
    }

    const readdy::scalar totalEnergy = energy03 + energy05 + energy13 + energy15 + energy23 + energy25;
    if(kernel->singlePrecision()) {
        EXPECT_FLOAT_EQ(stateModel.energy(), totalEnergy);
    } else {
        EXPECT_DOUBLE_EQ(stateModel.energy(), totalEnergy);
    }
}

TEST_P(TestStateModel, CalculateForcesNoForces) {
    m::Context &ctx = kernel->context();
    auto &stateModel = kernel->stateModel();
    auto calculateForces = kernel->actions().calculateForces();
    // several particles without potentials -> forces must all be zero
    auto obs = kernel->observe().forces(1);
    auto conn = kernel->connectObservable(obs.get());
    ctx.particleTypes().add("A", 1.0);
    ctx.particleTypes().add("B", 1.0);
    ctx.boxSize() = {{4., 4., 4.}};
    ctx.periodicBoundaryConditions() = {{false, false, false}};
    ctx.potentials().addBox("A", .00000001, {-1.9, -1.9, -1.9}, {3.85, 3.85, 3.85});
    ctx.potentials().addBox("B", .00000001, {-1.9, -1.9, -1.9}, {3.85, 3.85, 3.85});
    auto typeIdA = ctx.particleTypes().idOf("A");
    auto typeIdB = ctx.particleTypes().idOf("B");
    auto particlesA = std::vector<m::Particle> {
            m::Particle(0, 0, 0, typeIdA), m::Particle(0, 0.8, 0, typeIdA), m::Particle(0.2, 0, -0.2, typeIdA)
    };
    auto particlesB = std::vector<m::Particle> {
            m::Particle(0, 0, 0.1, typeIdB), m::Particle(-1.8, 0, 1.8, typeIdB), m::Particle(0.5, 0, 0, typeIdB)
    };
    stateModel.addParticles(particlesA);
    stateModel.addParticles(particlesB);
    stateModel.initializeNeighborList(0.);
    stateModel.updateNeighborList();
    calculateForces->perform();
    // check results
    obs->evaluate();
    for (auto &&force : obs->getResult()) {
        EXPECT_VEC3_EQ(force, readdy::Vec3(0, 0, 0));
    }
}

INSTANTIATE_TEST_CASE_P(TestStateModel, TestStateModel,
                        ::testing::ValuesIn(readdy::testing::getKernelsToTest()));
}
