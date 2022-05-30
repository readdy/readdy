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
 * @file TestStateModel.cpp
 * @brief Test the methods that manipulate time-dependent simulation data for kernels.
 * @author clonker
 * @author chrisfroe
 * @date 25.08.16
 * @todo check force calculation through periodic boundary
 */

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <readdy/testing/Utils.h>
#include <readdy/testing/KernelTest.h>

namespace m = readdy::model;

using namespace readdytesting::kernel;

TEMPLATE_TEST_CASE("Test state model", "[state-model]", SingleCPU, CPU) {
    auto kernel = create<TestType>();
    SECTION("Calculate force between two particles") {
        m::Context &ctx = kernel->context();
        auto &stateModel = kernel->stateModel();

        auto obs = kernel->observe().forces(1);
        auto conn = kernel->connectObservable(obs.get());
        // two A particles with radius 1. -> cutoff 2, distance 1.8 -> r-r_0 = 0.2 -> force = 0.2
        ctx.particleTypes().add("A", 1.0);
        ctx.boxSize() = {{4., 4., 4.}};
        ctx.periodicBoundaryConditions() = {{false, false, false}};

        ctx.potentials().addBox("A", .001, {-1.9, -1.9, -1.9}, {3.8, 3.8, 3.8});

        ctx.potentials().addHarmonicRepulsion("A", "A", 1.0, 2.0);
        ctx.potentials().addBox("A", .01, {-1.9, -1.9, -1.9}, {3.8, 3.8, 3.8});

        kernel->initialize();
        auto typeIdA = ctx.particleTypes().idOf("A");
        auto twoParticles = std::vector<m::Particle> {m::Particle(0., 0., 0., typeIdA), m::Particle(0., 0., 1.8, typeIdA)};

        stateModel.addParticles(twoParticles);
        stateModel.initializeNeighborList(ctx.calculateMaxCutoff());
        stateModel.updateNeighborList();

        auto calculateForces = kernel->actions().calculateForces();
        calculateForces->perform();
        calculateForces->perform(); // calculating twice should yield the same result. force and energy must not accumulate
        // check results
        obs->evaluate();
        auto forcesIt = obs->getResult().begin();
        if(kernel->doublePrecision()) {
            readdy::testing::vec3eq(*forcesIt, readdy::Vec3(0, 0, -0.2));
        } else {
            readdy::testing::vec3eq(*forcesIt, readdy::Vec3(0, 0, -0.2));
        }
        ++forcesIt;
        if(kernel->doublePrecision()) {
            readdy::testing::vec3eq(*forcesIt, readdy::Vec3(0, 0, 0.2));
        } else {
            readdy::testing::vec3eq(*forcesIt, readdy::Vec3(0, 0, 0.2));
        }
        REQUIRE(stateModel.energy() == Catch::Approx(0.02));
    }

    SECTION("Calculate repulsion forces") {
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
        std::vector<readdy::ParticleId> ids;
        {
            std::for_each(particlesA.begin(), particlesA.end(), [&ids](const m::Particle& p) { ids.push_back(p.id());});
            std::for_each(particlesB.begin(), particlesB.end(), [&ids](const m::Particle& p) { ids.push_back(p.id());});
        }
        stateModel.addParticles(particlesA);
        stateModel.addParticles(particlesB);
        kernel->initialize();
        stateModel.initializeNeighborList(ctx.calculateMaxCutoff());
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
            if(particle.id() == ids.at(0)) {
                readdy::testing::vec3eq(forces.at(idx), force03 + force05);
            } else if(particle.id() == ids.at(1)) {
                readdy::testing::vec3eq(forces.at(idx), force13 + force15);
            } else if(particle.id() == ids.at(2)) {
                readdy::testing::vec3eq(forces.at(idx), force23 + force25);
            } else if(particle.id() == ids.at(3)) {
                readdy::testing::vec3eq(forces.at(idx), (-1. * force03) - force13 - force23);
            } else if(particle.id() == ids.at(4)) {
                readdy::testing::vec3eq(forces.at(idx), readdy::Vec3(0, 0, 0));
            } else if(particle.id() == ids.at(5)) {
                readdy::testing::vec3eq(forces.at(idx), (-1. * force05) - force15 - force25);
            } else {
                readdy::log::error("Got an unexpected particle id: {}", particle.id());
                FAIL();
            }
            ++idx;
        }

        const readdy::scalar totalEnergy = energy03 + energy05 + energy13 + energy15 + energy23 + energy25;
        REQUIRE(stateModel.energy() == Catch::Approx(totalEnergy));
    }

    SECTION("Calculate no forces") {
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
        calculateForces->perform();
        // check results
        obs->evaluate();
        for (auto &&force : obs->getResult()) {
            readdy::testing::vec3eq(force, readdy::Vec3(0, 0, 0));
        }
    }
}
