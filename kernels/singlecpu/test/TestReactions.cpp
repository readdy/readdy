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
 * << detailed description >>
 *
 * @file TestReactions.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 22.06.16
 */

#include <gtest/gtest.h>
#include <readdy/model/Kernel.h>
#include <readdy/plugin/KernelProvider.h>
#include <readdy/model/actions/Actions.h>

TEST(SingleCPUTestReactions, TestDecay) {
    using particle_t = readdy::model::Particle;
    auto kernel = readdy::plugin::KernelProvider::getInstance().create("SingleCPU");
    kernel->context().boxSize() = {{10, 10, 10}};
    kernel->context().particleTypes().add("X", .25);
    kernel->context().reactions().addDecay("X decay", "X", 1e16);
    kernel->context().reactions().addFission("X fission", "X", "X", "X", .5, .3);

    readdy::scalar timeStep = 1.0;
    auto &&integrator = kernel->actions().eulerBDIntegrator(timeStep);
    auto &&forces = kernel->actions().calculateForces();
    using update_nl = readdy::model::actions::UpdateNeighborList;
    auto &&initNeighborList = kernel->actions().updateNeighborList(update_nl::Operation::init, 0);
    auto &&neighborList = kernel->actions().updateNeighborList(update_nl::Operation::update, 0);
    auto &&reactions = kernel->actions().uncontrolledApproximation(timeStep);

    auto pp_obs = kernel->observe().positions(1);
    pp_obs->callback() = [](const readdy::model::observables::Positions::result_type &t) {
        readdy::log::trace("got n particles={}", t.size());
    };
    auto connection = kernel->connectObservable(pp_obs.get());

    const int n_particles = 200;
    const auto typeId = kernel->context().particleTypes().idOf("X");
    std::vector<readdy::model::Particle> particlesToBeginWith{n_particles, {0, 0, 0, typeId}};
    kernel->stateModel().addParticles(particlesToBeginWith);
    initNeighborList->perform();
    for (size_t t = 0; t < 20; t++) {

        forces->perform();
        integrator->perform();
        neighborList->perform();
        reactions->perform();

        kernel->evaluateObservables(t);

    }

    EXPECT_EQ(0, kernel->stateModel().getParticlePositions().size());

    connection.disconnect();
}

/**
 * Setting:
 *  - Five particle types A,B,C,D,E.
 *  - They are instantiated such that there is one A particle, one B and one C particle.
 *  - B and C particle are within the reaction radius of their assigned reaction.
 * Reactions (all with rate 1, dt = 1):
 *  - A has no interaction with the other particles and dies after one time step (death rate 1)
 *  - B + C -> E
 *  - B + D -> A
 *  - E -> A
 *  - C -> D
 * Expected:
 *  - After one time step, the A particle dies and either C -> D or B + C -> E
 *  - After two time steps:
 *      - If previously C -> D, one is left with one B and one D particle,
 *        which will react to one A particle
 *      - If previously B + C -> E, one is left with one E particle,
 *        which will convert to one A particle
 *  - After three time steps, one is left with one A particle which then will
 *    decay within the next timestep, leaving no particles.
 * Assert:
 *   - t = 0: n_particles == 3 with 1x A, 1x B, 1x C
 *   - t = 1: n_particles == 2 || n_particles == 1 with (1x B, 1x D) || 1x E
 *   - t = 2: n_particles == 1 with 1x A
 *   - t > 2: n_particles == 0
 */
TEST(SingleCPUTestReactions, TestMultipleReactionTypes) {
    using particle_t = readdy::model::Particle;
    auto kernel = readdy::plugin::KernelProvider::getInstance().create("SingleCPU");
    kernel->context().boxSize() = {{10, 10, 10}};

    kernel->context().particleTypes().add("A", .25);
    kernel->context().particleTypes().add("B", .25);
    kernel->context().particleTypes().add("C", .25);
    kernel->context().particleTypes().add("D", .25);
    kernel->context().particleTypes().add("E", .25);

    kernel->context().reactions().addDecay("A decay", "A", 1e16);
    kernel->context().reactions().addFusion("B+C->E", "B", "C", "E", 1e16, 17);
    kernel->context().reactions().addFusion("B+D->A", "B", "D", "A", 1e16, 17);
    kernel->context().reactions().addConversion("E->A", "E", "A", 1e16);
    kernel->context().reactions().addConversion("C->D", "C", "D", 1e16);

    auto &&integrator = kernel->actions().eulerBDIntegrator(1);
    auto &&forces = kernel->actions().calculateForces();
    auto &&neighborList = kernel->actions().updateNeighborList();
    auto &&reactions = kernel->actions().uncontrolledApproximation(1);

    const auto typeId_A = kernel->context().particleTypes().idOf("A");
    const auto typeId_B = kernel->context().particleTypes().idOf("B");
    const auto typeId_C = kernel->context().particleTypes().idOf("C");
    const auto typeId_D = kernel->context().particleTypes().idOf("D");
    const auto typeId_E = kernel->context().particleTypes().idOf("E");

    kernel->stateModel().addParticle({4, 4, 4, typeId_A});
    kernel->stateModel().addParticle({-2, 0, 0, typeId_B});
    kernel->stateModel().addParticle({2, 0, 0, typeId_C});

    auto pred_contains_A = [=](const readdy::model::Particle &p) { return p.type() == typeId_A; };
    auto pred_contains_B = [=](const readdy::model::Particle &p) { return p.type() == typeId_B; };
    auto pred_contains_C = [=](const readdy::model::Particle &p) { return p.type() == typeId_C; };
    auto pred_contains_D = [=](const readdy::model::Particle &p) { return p.type() == typeId_D; };
    auto pred_contains_E = [=](const readdy::model::Particle &p) { return p.type() == typeId_E; };

    for (unsigned int t = 0; t < 4; t++) {

        const auto particles = kernel->stateModel().getParticles();

        bool containsA = std::find_if(particles.begin(), particles.end(), pred_contains_A) != particles.end();
        bool containsB = std::find_if(particles.begin(), particles.end(), pred_contains_B) != particles.end();
        bool containsC = std::find_if(particles.begin(), particles.end(), pred_contains_C) != particles.end();
        bool containsD = std::find_if(particles.begin(), particles.end(), pred_contains_D) != particles.end();
        bool containsE = std::find_if(particles.begin(), particles.end(), pred_contains_E) != particles.end();

        switch (t) {
            case 0: {
                EXPECT_EQ(3, particles.size());
                EXPECT_TRUE(containsA);
                EXPECT_TRUE(containsB);
                EXPECT_TRUE(containsC);
                break;
            }
            case 1: {
                EXPECT_TRUE(particles.size() == 2 || particles.size() == 1);
                if (particles.size() == 2) {
                    readdy::log::debug("------> conversion happened");
                    EXPECT_TRUE(containsB);
                    EXPECT_TRUE(containsD);
                } else {
                    readdy::log::debug("------> fusion happened");
                    EXPECT_TRUE(containsE);
                }
                break;
            }
            case 2: {
                EXPECT_EQ(1, particles.size());
                EXPECT_TRUE(containsA);
                break;
            }
            case 3: {
                EXPECT_EQ(0, particles.size());
                break;
            }
            default: {
                FAIL();
            }
        }

        // propagate
        neighborList->perform();
        forces->perform();
        integrator->perform();

        neighborList->perform();
        reactions->perform();
        kernel->evaluateObservables(t);
    }
}
