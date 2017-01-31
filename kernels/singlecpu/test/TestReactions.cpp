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

TEST(SingleCPUTestReactions, CheckInOutTypesAndPositions) {
    using fusion_t = readdy::model::reactions::Fusion;
    using fission_t = readdy::model::reactions::Fission;
    using enzymatic_t = readdy::model::reactions::Enzymatic;
    using conversion_t = readdy::model::reactions::Conversion;
    using death_t = readdy::model::reactions::Decay;
    using particle_t = readdy::model::Particle;
    auto kernel = readdy::plugin::KernelProvider::getInstance().create("SingleCPU");
    kernel->getKernelContext().setPeriodicBoundary(false, false, false);
    kernel->getKernelContext().setBoxSize(100, 100, 100);
    const auto diff = kernel->getKernelContext().getShortestDifferenceFun();
    kernel->getKernelContext().registerParticleType("A", .1, 1.); // type id 0
    kernel->getKernelContext().registerParticleType("B", .1, 1.); // type id 1
    kernel->getKernelContext().registerParticleType("C", .1, 1.); // type id 2

    // test conversion
    {
        auto conversion = kernel->getReactionFactory().createReaction<conversion_t>("A->B", 0, 1, 1);
        particle_t p_A{0, 0, 0, 0};
        particle_t p_out{5, 5, 5, 1};
        conversion->perform(p_A, p_A, p_out, p_out);
        EXPECT_EQ(p_out.getType(), conversion->getTypeTo());
        EXPECT_EQ(p_out.getPos(), p_A.getPos());
    }

    // test fusion
    {
        double eductDistance = .4;
        double weight1 = .3, weight2 = .7;
        auto fusion = kernel->getReactionFactory().createReaction<fusion_t>("A+B->C", 0, 1, 2, 1, eductDistance,
                                                                            weight1, weight2);
        particle_t p_out1{50, 50, 50, 70};
        particle_t p_out2{50, 50, 50, 70};
        particle_t p_A{1, 0, 0, 0};
        particle_t p_B{-1, 0, 0, 1};
        fusion->perform(p_A, p_B, p_out1, p_out1);
        fusion->perform(p_B, p_A, p_out2, p_out2);

        EXPECT_EQ(p_out1.getPos(), p_out2.getPos());
        EXPECT_EQ(p_out1.getType(), p_out2.getType());
        EXPECT_EQ(p_out1.getType(), fusion->getTo());

        EXPECT_EQ(readdy::model::Vec3(.4, 0, 0), p_out1.getPos());
    }

    // fission
    {
        double productDistance = .4;
        double weight1 = .3, weight2 = .7;
        auto fission = kernel->getReactionFactory().createReaction<fission_t>("C->A+B", 2, 0, 1, 1, productDistance,
                                                                              weight1, weight2);
        particle_t p_C{0, 0, 0, 2};
        particle_t p_out1{50, 50, 50, 70};
        particle_t p_out2{50, 50, 50, 70};
        fission->perform(p_C, p_C, p_out1, p_out2);

        EXPECT_EQ(p_out1.getType(), fission->getTo1());
        EXPECT_EQ(p_out2.getType(), fission->getTo2());
        auto p_12 = diff(p_out1.getPos(), p_out2.getPos());
        auto p_12_nondirect = p_out2.getPos() - p_out1.getPos();
        EXPECT_EQ(p_12_nondirect, p_12);
        auto distance = std::sqrt(p_12 * p_12);
        EXPECT_DOUBLE_EQ(productDistance, distance);
    }

    // enzymatic
    {
        auto enzymatic = kernel->getReactionFactory().createReaction<enzymatic_t>("A+C->B+C", 2, 0, 1, 1, .5);
        particle_t p_A{0, 0, 0, 0};
        particle_t p_C{5, 5, 5, 2};
        {
            particle_t p_out1{50, 50, 50, 70};
            particle_t p_out2{50, 50, 50, 70};
            enzymatic->perform(p_A, p_C, p_out1, p_out2);
            if (p_out1.getType() == enzymatic->getCatalyst()) {
                EXPECT_EQ(enzymatic->getCatalyst(), p_out1.getType());
                EXPECT_EQ(enzymatic->getTo(), p_out2.getType());
                EXPECT_EQ(p_C.getPos(), p_out1.getPos());
                EXPECT_EQ(p_A.getPos(), p_out2.getPos());
            } else {
                EXPECT_EQ(enzymatic->getCatalyst(), p_out2.getType());
                EXPECT_EQ(enzymatic->getTo(), p_out1.getType());
                EXPECT_EQ(p_C.getPos(), p_out2.getPos());
                EXPECT_EQ(p_A.getPos(), p_out1.getPos());
            }
        }
        {
            particle_t p_out1{50, 50, 50, 70};
            particle_t p_out2{50, 50, 50, 70};
            enzymatic->perform(p_C, p_A, p_out1, p_out2);
            if (p_out1.getType() == enzymatic->getCatalyst()) {
                EXPECT_EQ(enzymatic->getCatalyst(), p_out1.getType());
                EXPECT_EQ(enzymatic->getTo(), p_out2.getType());
                EXPECT_EQ(p_C.getPos(), p_out1.getPos());
                EXPECT_EQ(p_A.getPos(), p_out2.getPos());
            } else {
                EXPECT_EQ(enzymatic->getCatalyst(), p_out2.getType());
                EXPECT_EQ(enzymatic->getTo(), p_out1.getType());
                EXPECT_EQ(p_C.getPos(), p_out2.getPos());
                EXPECT_EQ(p_A.getPos(), p_out1.getPos());
            }
        }
    }

}

TEST(SingleCPUTestReactions, TestDecay) {
    using fusion_t = readdy::model::reactions::Fusion;
    using fission_t = readdy::model::reactions::Fission;
    using enzymatic_t = readdy::model::reactions::Enzymatic;
    using conversion_t = readdy::model::reactions::Conversion;
    using death_t = readdy::model::reactions::Decay;
    using particle_t = readdy::model::Particle;
    auto kernel = readdy::plugin::KernelProvider::getInstance().create("SingleCPU");
    kernel->getKernelContext().setBoxSize(10, 10, 10);
    kernel->getKernelContext().registerParticleType("X", .25, 1.);
    kernel->registerReaction<death_t>("X decay", "X", 1);
    kernel->registerReaction<fission_t>("X fission", "X", "X", "X", .5, .3);

    double timeStep = 1.0;
    auto &&integrator = kernel->createAction<readdy::model::actions::EulerBDIntegrator>(timeStep);
    auto &&forces = kernel->createAction<readdy::model::actions::CalculateForces>();
    using update_nl = readdy::model::actions::UpdateNeighborList;
    auto &&neighborList = kernel->createAction<readdy::model::actions::UpdateNeighborList>(update_nl::Operation::create, -1);
    auto &&reactions = kernel->createAction<readdy::model::actions::reactions::UncontrolledApproximation>(timeStep);

    auto pp_obs = kernel->createObservable<readdy::model::observables::Positions>(1);
    pp_obs->setCallback([](const readdy::model::observables::Positions::result_t &t) {
        readdy::log::console()->trace("got n particles={}", t.size());
    });
    auto connection = kernel->connectObservable(pp_obs.get());

    const int n_particles = 200;
    const auto typeId = kernel->getKernelContext().getParticleTypeID("X");
    std::vector<readdy::model::Particle> particlesToBeginWith{n_particles, {0, 0, 0, typeId}};
    kernel->getKernelStateModel().addParticles(particlesToBeginWith);
    kernel->getKernelContext().configure();
    neighborList->perform();
    for (size_t t = 0; t < 20; t++) {

        forces->perform();
        integrator->perform();
        neighborList->perform();
        reactions->perform();

        kernel->evaluateObservables(t);

    }

    EXPECT_EQ(0, kernel->getKernelStateModel().getParticlePositions().size());

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
    using fusion_t = readdy::model::reactions::Fusion;
    using fission_t = readdy::model::reactions::Fission;
    using enzymatic_t = readdy::model::reactions::Enzymatic;
    using conversion_t = readdy::model::reactions::Conversion;
    using death_t = readdy::model::reactions::Decay;
    using particle_t = readdy::model::Particle;
    auto kernel = readdy::plugin::KernelProvider::getInstance().create("SingleCPU");
    kernel->getKernelContext().setBoxSize(10, 10, 10);

    kernel->getKernelContext().registerParticleType("A", .25, 1.);
    kernel->getKernelContext().registerParticleType("B", .25, 1.);
    kernel->getKernelContext().registerParticleType("C", .25, 1.);
    kernel->getKernelContext().registerParticleType("D", .25, 1.);
    kernel->getKernelContext().registerParticleType("E", .25, 1.);

    kernel->registerReaction<death_t>("A decay", "A", 1);
    kernel->registerReaction<fusion_t>("B+C->E", "B", "C", "E", 1, 17);
    kernel->registerReaction<fusion_t>("B+D->A", "B", "D", "A", 1, 17);
    kernel->registerReaction<conversion_t>("E->A", "E", "A", 1);
    kernel->registerReaction<conversion_t>("C->D", "C", "D", 1);

    auto &&integrator = kernel->createAction<readdy::model::actions::EulerBDIntegrator>(1);
    auto &&forces = kernel->createAction<readdy::model::actions::CalculateForces>();
    auto &&neighborList = kernel->createAction<readdy::model::actions::UpdateNeighborList>();
    auto &&reactions = kernel->createAction<readdy::model::actions::reactions::UncontrolledApproximation>(1);

    const auto typeId_A = kernel->getKernelContext().getParticleTypeID("A");
    const auto typeId_B = kernel->getKernelContext().getParticleTypeID("B");
    const auto typeId_C = kernel->getKernelContext().getParticleTypeID("C");
    const auto typeId_D = kernel->getKernelContext().getParticleTypeID("D");
    const auto typeId_E = kernel->getKernelContext().getParticleTypeID("E");

    kernel->getKernelStateModel().addParticle({4, 4, 4, typeId_A});
    kernel->getKernelStateModel().addParticle({-2, 0, 0, typeId_B});
    kernel->getKernelStateModel().addParticle({2, 0, 0, typeId_C});

    auto pred_contains_A = [=](const readdy::model::Particle &p) { return p.getType() == typeId_A; };
    auto pred_contains_B = [=](const readdy::model::Particle &p) { return p.getType() == typeId_B; };
    auto pred_contains_C = [=](const readdy::model::Particle &p) { return p.getType() == typeId_C; };
    auto pred_contains_D = [=](const readdy::model::Particle &p) { return p.getType() == typeId_D; };
    auto pred_contains_E = [=](const readdy::model::Particle &p) { return p.getType() == typeId_E; };

    kernel->getKernelContext().configure();

    for (unsigned int t = 0; t < 4; t++) {

        const auto particles = kernel->getKernelStateModel().getParticles();

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
                    readdy::log::console()->debug("------> conversion happened");
                    EXPECT_TRUE(containsB);
                    EXPECT_TRUE(containsD);
                } else {
                    readdy::log::console()->debug("------> fusion happened");
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
