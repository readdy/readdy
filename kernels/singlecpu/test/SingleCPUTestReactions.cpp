/**
 * << detailed description >>
 *
 * @file TestReactions.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 22.06.16
 */

#include <gtest/gtest.h>
#include <readdy/kernel/singlecpu/reactions/SingleCPUReactions.h>
#include <readdy/model/Kernel.h>
#include <boost/algorithm/string.hpp>
#include <readdy/plugin/KernelProvider.h>
#include <readdy/model/programs/Programs.h>

TEST(SingleCPUTestReactions, CheckInOutTypesAndPositions) {
    using fusion_t = readdy::model::reactions::Fusion;
    using fission_t = readdy::model::reactions::Fission;
    using enzymatic_t = readdy::model::reactions::Enzymatic;
    using conversion_t = readdy::model::reactions::Conversion;
    using death_t = readdy::model::reactions::Death;
    using particle_t = readdy::model::Particle;
    auto kernel = readdy::plugin::KernelProvider::getInstance().create("SingleCPU");
    kernel->getKernelContext().setPeriodicBoundary(false, false, false);
    kernel->getKernelContext().setBoxSize(100, 100, 100);
    const auto diff = kernel->getKernelContext().getShortestDifferenceFun();
    kernel->getKernelContext().setDiffusionConstant("A", .1); // type id 0
    kernel->getKernelContext().setDiffusionConstant("B", .1); // type id 1
    kernel->getKernelContext().setDiffusionConstant("C", .1); // type id 2

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
        auto fusion = kernel->getReactionFactory().createReaction<fusion_t>("A+B->C", 0, 1, 2, 1, eductDistance, weight1, weight2);
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
        auto fission = kernel->getReactionFactory().createReaction<fission_t>("C->A+B", 2, 0, 1, productDistance, 1, weight1, weight2);
        particle_t p_C{0, 0, 0, 2};
        particle_t p_out1{50, 50, 50, 70};
        particle_t p_out2{50, 50, 50, 70};
        fission->perform(p_C, p_C, p_out1, p_out2);

        EXPECT_EQ(p_out1.getType(), fission->getTo1());
        EXPECT_EQ(p_out2.getType(), fission->getTo2());
        auto p_12 = diff(p_out1.getPos(), p_out2.getPos());
        auto p_12_nondirect = p_out2.getPos() - p_out1.getPos();
        EXPECT_EQ(p_12_nondirect, p_12);
        auto distance = sqrt(p_12 * p_12);
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
    using death_t = readdy::model::reactions::Death;
    using particle_t = readdy::model::Particle;
    auto kernel = readdy::plugin::KernelProvider::getInstance().create("SingleCPU");
    kernel->getKernelContext().setBoxSize(10, 10, 10);
    kernel->getKernelContext().setTimeStep(1);
    kernel->getKernelContext().setDiffusionConstant("X", .25);
    kernel->getKernelContext().registerDeathReaction("X decay", "X", .5);
    kernel->getKernelContext().registerFissionReaction("X fission", "X", "X", "X", .15, .00);

    auto &&integrator = kernel->createProgram<readdy::model::programs::EulerBDIntegrator>();
    auto &&forces = kernel->createProgram<readdy::model::programs::CalculateForces>();
    auto &&neighborList = kernel->createProgram<readdy::model::programs::UpdateNeighborList>();
    auto &&reactionsProgram = kernel->createProgram<readdy::model::programs::reactions::UncontrolledApproximation>();

    auto pp_obs = kernel->createObservable<readdy::model::ParticlePositionObservable>(1);
    auto connection = kernel->connectObservable(pp_obs.get());

    const int n_particles = 200;
    const unsigned int typeId = kernel->getKernelContext().getParticleTypeID("X");
    std::vector<readdy::model::Particle> particlesToBeginWith{n_particles, {0, 0, 0, typeId}};
    kernel->getKernelStateModel().addParticles(particlesToBeginWith);

    neighborList->execute();
    for (size_t t = 0; t < 20; t++) {

        forces->execute();
        integrator->execute();
        neighborList->execute();
        reactionsProgram->execute();

        pp_obs->evaluate();

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
    using death_t = readdy::model::reactions::Death;
    using particle_t = readdy::model::Particle;
    auto kernel = readdy::plugin::KernelProvider::getInstance().create("SingleCPU");
    kernel->getKernelContext().setBoxSize(10, 10, 10);
    kernel->getKernelContext().setTimeStep(1);

    kernel->getKernelContext().setDiffusionConstant("A", .25);
    kernel->getKernelContext().setDiffusionConstant("B", .25);
    kernel->getKernelContext().setDiffusionConstant("C", .25);
    kernel->getKernelContext().setDiffusionConstant("D", .25);
    kernel->getKernelContext().setDiffusionConstant("E", .25);

    kernel->getKernelContext().registerDeathReaction("A decay", "A", 1);
    kernel->getKernelContext().registerFusionReaction("B+C->E", "B", "C", "E", 1, 13);
    kernel->getKernelContext().registerFusionReaction("B+D->A", "B", "D", "A", 1, 13);
    kernel->getKernelContext().registerConversionReaction("E->A", "E", "A", 1);
    kernel->getKernelContext().registerConversionReaction("C->D", "C", "D", 1);

    auto &&integrator = kernel->createProgram<readdy::model::programs::EulerBDIntegrator>();
    auto &&forces = kernel->createProgram<readdy::model::programs::CalculateForces>();
    auto &&neighborList = kernel->createProgram<readdy::model::programs::UpdateNeighborList>();
    auto &&reactionsProgram = kernel->createProgram<readdy::model::programs::reactions::UncontrolledApproximation>();

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
                    BOOST_LOG_TRIVIAL(debug) << "------> conversion happened";
                    EXPECT_TRUE(containsB);
                    EXPECT_TRUE(containsD);
                } else {
                    BOOST_LOG_TRIVIAL(debug) << "------> fusion happened";
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
        neighborList->execute();
        forces->execute();
        integrator->execute();

        neighborList->execute();
        reactionsProgram->execute();
    }
}
