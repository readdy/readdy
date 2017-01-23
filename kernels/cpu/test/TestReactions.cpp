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
 * @date 01.09.16
 */

#include <gtest/gtest.h>
#include <readdy/model/Kernel.h>
#include <readdy/plugin/KernelProvider.h>
#include <readdy/kernel/cpu/CPUKernel.h>
#include <readdy/kernel/cpu/actions/reactions/ReactionUtils.h>
#include <readdy/kernel/cpu/actions/reactions/CPUGillespieParallel.h>
#include <readdy/kernel/cpu/actions/reactions/NextSubvolumesReactionScheduler.h>

namespace reac = readdy::kernel::cpu::actions::reactions;

struct fix_n_threads {
    fix_n_threads(readdy::kernel::cpu::CPUKernel *const kernel, unsigned int n)
            : oldValue(static_cast<unsigned int>(kernel->getNThreads())), kernel(kernel) {
        kernel->setNThreads(n);
    }

    ~fix_n_threads() {
        kernel->setNThreads(oldValue);
    }

private:
    const unsigned int oldValue;
    readdy::kernel::cpu::CPUKernel *const kernel;
};

TEST(CPUTestReactions, CheckInOutTypesAndPositions) {
    using fusion_t = readdy::model::reactions::Fusion;
    using fission_t = readdy::model::reactions::Fission;
    using enzymatic_t = readdy::model::reactions::Enzymatic;
    using conversion_t = readdy::model::reactions::Conversion;
    using death_t = readdy::model::reactions::Decay;
    using particle_t = readdy::model::Particle;
    using data_t = readdy::kernel::cpu::model::CPUParticleData;
    auto kernel = readdy::plugin::KernelProvider::getInstance().create("CPU");
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

        data_t data {&kernel->getKernelContext()};
        data.addParticles({p_A});

        data_t::entries_update_t newParticles{};
        std::vector<data_t::index_t> decayedEntries {};

        reac::performReaction(data, 0, 0, newParticles, decayedEntries, conversion.get());

        EXPECT_EQ(data.entry_at(0).type, conversion->getTypeTo());
        EXPECT_EQ(data.pos(0), readdy::model::Vec3(0,0,0));
    }

    // test fusion
    {
        data_t data {&kernel->getKernelContext()};

        data_t::entries_update_t newParticles{};
        std::vector<data_t::index_t> decayedEntries {};

        double eductDistance = .4;
        double weight1 = .3, weight2 = .7;
        auto fusion = kernel->getReactionFactory().createReaction<fusion_t>("A+B->C", 0, 1, 2, 1, eductDistance,
                                                                            weight1, weight2);
        particle_t p_A{1, 0, 0, 0};
        particle_t p_B{-1, 0, 0, 1};
        data.addParticles({p_A, p_B});

        reac::performReaction(data, 0, 1, newParticles, decayedEntries, fusion.get());

        EXPECT_EQ(data.entry_at(0).type, fusion->getTo());
        EXPECT_EQ(decayedEntries.size(), 1);
        EXPECT_EQ(decayedEntries.at(0), 1);

        EXPECT_EQ(readdy::model::Vec3(.4, 0, 0), data.pos(0));
    }

    // fission
    {
        data_t data {&kernel->getKernelContext()};

        data_t::entries_update_t newParticles{};
        std::vector<data_t::index_t> decayedEntries {};

        double productDistance = .4;
        double weight1 = .3, weight2 = .7;
        auto fission = kernel->getReactionFactory().createReaction<fission_t>("C->A+B", 2, 0, 1, 1, productDistance,
                                                                              weight1, weight2);
        particle_t p_C{0, 0, 0, 2};
        data.addParticle(p_C);

        reac::performReaction(data, 0, 0, newParticles, decayedEntries, fission.get());
        data.update(std::make_pair(std::move(newParticles), std::move(decayedEntries)));

        EXPECT_EQ(data.entry_at(0).type, fission->getTo1());
        EXPECT_EQ(data.entry_at(1).type, fission->getTo2());
        auto p_12 = diff(data.pos(0), data.pos(1));
        auto p_12_nondirect = data.pos(1) - data.pos(0);
        EXPECT_EQ(p_12_nondirect, p_12);
        auto distance = std::sqrt(p_12 * p_12);
        EXPECT_DOUBLE_EQ(productDistance, distance);
    }

    // enzymatic 1
    {
        data_t data{&kernel->getKernelContext()};

        data_t::entries_update_t newParticles{};
        std::vector<data_t::index_t> decayedEntries{};

        auto enzymatic = kernel->getReactionFactory().createReaction<enzymatic_t>("A+C->B+C", 2, 0, 1, 1, .5);
        particle_t p_A{0, 0, 0, 0};
        particle_t p_C{5, 5, 5, 2};
        data.addParticles({p_A, p_C});
        reac::performReaction(data, 0, 1, newParticles, decayedEntries, enzymatic.get());
        data.update(std::make_pair(std::move(newParticles), std::move(decayedEntries)));
        {
            const auto &e1 = data.entry_at(0);
            const auto &e2 = data.entry_at(1);
            if (e1.type == enzymatic->getCatalyst()) {
                EXPECT_EQ(enzymatic->getCatalyst(), e1.type);
                EXPECT_EQ(enzymatic->getTo(), e2.type);
                EXPECT_EQ(p_C.getPos(), e1.position());
                EXPECT_EQ(p_A.getPos(), e2.position());
            } else {
                EXPECT_EQ(enzymatic->getCatalyst(), e2.type);
                EXPECT_EQ(enzymatic->getTo(), e1.type);
                EXPECT_EQ(p_C.getPos(), e2.position());
                EXPECT_EQ(p_A.getPos(), e1.position());
            }
        }
    }
    // enzymatic 2
    {
        data_t data{&kernel->getKernelContext()};

        data_t::entries_update_t newParticles{};
        std::vector<data_t::index_t> decayedEntries{};

        auto enzymatic = kernel->getReactionFactory().createReaction<enzymatic_t>("A+C->B+C", 2, 0, 1, 1, .5);
        particle_t p_A{0, 0, 0, 0};
        particle_t p_C{5, 5, 5, 2};
        data.addParticles({p_C, p_A});
        reac::performReaction(data, 0, 1, newParticles, decayedEntries, enzymatic.get());
        data.update(std::make_pair(std::move(newParticles), std::move(decayedEntries)));
        {
            const auto &e1 = data.entry_at(0);
            const auto &e2 = data.entry_at(1);
            if (e1.type == enzymatic->getCatalyst()) {
                EXPECT_EQ(enzymatic->getCatalyst(), e1.type);
                EXPECT_EQ(enzymatic->getTo(), e2.type);
                EXPECT_EQ(p_C.getPos(), e1.position());
                EXPECT_EQ(p_A.getPos(), e2.position());
            } else {
                EXPECT_EQ(enzymatic->getCatalyst(), e2.type);
                EXPECT_EQ(enzymatic->getTo(), e1.type);
                EXPECT_EQ(p_C.getPos(), e2.position());
                EXPECT_EQ(p_A.getPos(), e1.position());
            }
        }
    }

}

TEST(CPUTestReactions, TestDecay) {
    using fusion_t = readdy::model::reactions::Fusion;
    using fission_t = readdy::model::reactions::Fission;
    using enzymatic_t = readdy::model::reactions::Enzymatic;
    using conversion_t = readdy::model::reactions::Conversion;
    using death_t = readdy::model::reactions::Decay;
    using particle_t = readdy::model::Particle;
    auto kernel = readdy::plugin::KernelProvider::getInstance().create("CPU");
    kernel->getKernelContext().setBoxSize(10, 10, 10);
    kernel->getKernelContext().setDiffusionConstant("X", .25);
    kernel->registerReaction<death_t>("X decay", "X", 1);
    kernel->registerReaction<fission_t>("X fission", "X", "X", "X", .5, .3);

    auto &&integrator = kernel->createAction<readdy::model::actions::EulerBDIntegrator>(1);
    auto &&forces = kernel->createAction<readdy::model::actions::CalculateForces>();
    auto &&neighborList = kernel->createAction<readdy::model::actions::UpdateNeighborList>();
    auto &&reactionsProgram = kernel->createAction<readdy::model::actions::reactions::GillespieParallel>(1);

    auto pp_obs = kernel->createObservable<readdy::model::observables::Positions>(1);
    pp_obs->setCallback([](const readdy::model::observables::Positions::result_t &t) {
        /* ignore */
    });
    auto connection = kernel->connectObservable(pp_obs.get());

    const int n_particles = 200;
    const unsigned int typeId = kernel->getKernelContext().getParticleTypeID("X");
    std::vector<readdy::model::Particle> particlesToBeginWith{n_particles, {0, 0, 0, typeId}};
    kernel->getKernelStateModel().addParticles(particlesToBeginWith);
    kernel->getKernelContext().configure();
    neighborList->perform();
    for (size_t t = 0; t < 20; t++) {

        forces->perform();
        integrator->perform();
        neighborList->perform();
        reactionsProgram->perform();

        kernel->evaluateObservables(t);

    }

    EXPECT_EQ(0, kernel->getKernelStateModel().getParticlePositions().size());

    connection.disconnect();
}

TEST(CPUTestReactions, TestGillespieParallel) {
    using particle_t = readdy::model::Particle;
    using vec_t = readdy::model::Vec3;
    using fusion_t = readdy::model::reactions::Fusion;
    using fission_t = readdy::model::reactions::Fission;
    using enzymatic_t = readdy::model::reactions::Enzymatic;
    using conversion_t = readdy::model::reactions::Conversion;
    using death_t = readdy::model::reactions::Decay;
    using particle_t = readdy::model::Particle;
    auto kernel = std::make_unique<readdy::kernel::cpu::CPUKernel>();
    kernel->getKernelContext().setBoxSize(10, 10, 30);
    kernel->getKernelContext().setPeriodicBoundary(true, true, false);

    kernel->getKernelContext().setDiffusionConstant("A", .25);
    kernel->getKernelContext().setDiffusionConstant("B", .25);
    kernel->getKernelContext().setDiffusionConstant("C", .25);
    double reactionRadius = 1.0;
    kernel->registerReaction<fusion_t>("annihilation", "A", "A", "A", 1.0, reactionRadius);
    kernel->registerReaction<fusion_t>("very unlikely", "A", "C", "A", std::numeric_limits<double>::min(),
                                       reactionRadius);
    kernel->registerReaction<fusion_t>("dummy reaction", "A", "B", "A", 0.0, reactionRadius);

    const auto typeA = kernel->getKernelContext().getParticleTypeID("A");
    const auto typeB = kernel->getKernelContext().getParticleTypeID("B");
    const auto typeC = kernel->getKernelContext().getParticleTypeID("C");

    // this particle goes right into the middle, i.e., into the halo region
    kernel->getKernelStateModel().addParticle({0, 0, 0, typeA});            // 0
    // these particles go left and right of this particle into the boxes as problematic ones
    kernel->getKernelStateModel().addParticle({0, 0, -.7, typeA});          // 1
    kernel->getKernelStateModel().addParticle({0, 0, .7, typeC});           // 2
    // these particles are far enough away from the halo region but still conflict (transitively) with the 1st layer
    kernel->getKernelStateModel().addParticle({0, 0, -1.6, typeC});         // 3
    kernel->getKernelStateModel().addParticle({0, 0, 1.6, typeA});          // 4
    kernel->getKernelStateModel().addParticle({0, 0, 1.7, typeA});          // 5
    // this particle are conflicting but should not appear as their reaction rate is 0
    kernel->getKernelStateModel().addParticle({0, 0, -1.7, typeB});         // 6
    // these particles are well inside the boxes and should not be considered problematic
    kernel->getKernelStateModel().addParticle({0, 0, -5, typeA});           // 7
    kernel->getKernelStateModel().addParticle({0, 0, -5.5, typeA});         // 8
    kernel->getKernelStateModel().addParticle({0, 0, 5, typeA});            // 9
    kernel->getKernelStateModel().addParticle({0, 0, 5.5, typeA});          // 10

    kernel->getKernelContext().configure();
    // a box width in z direction of 12 should divide into two boxes of 5x5x6
    {
        fix_n_threads n_threads{kernel.get(), 2};
        EXPECT_EQ(2, kernel->getNThreads());
        auto &&neighborList = kernel->createAction<readdy::model::actions::UpdateNeighborList>();
        auto &&reactionsProgram = readdy::util::static_unique_ptr_cast<readdy::kernel::cpu::actions::reactions::CPUGillespieParallel>(
                kernel->createAction<readdy::model::actions::reactions::GillespieParallel>(1)
        );
        neighborList->perform();
        reactionsProgram->perform();
        EXPECT_EQ(1.0, reactionsProgram->getMaxReactionRadius());
        EXPECT_EQ(15.0, reactionsProgram->getBoxWidth());
        EXPECT_EQ(2, reactionsProgram->getLongestAxis());
        EXPECT_TRUE(reactionsProgram->getOtherAxis1() == 0 || reactionsProgram->getOtherAxis1() == 1);
        if (reactionsProgram->getOtherAxis1() == 0) {
            EXPECT_EQ(1, reactionsProgram->getOtherAxis2());
        } else {
            EXPECT_EQ(0, reactionsProgram->getOtherAxis2());
        }
    }
    // we have two boxes, left and right, which can be projected into a line with (index, type):
    // Box left:  |---(8, A)---(7, A)--------(6, B)---(3, C)---(1, A)---(0, A)    |
    // Box right: |   (0, A)---(2, C)---(4, A)---(5, A)--------(9, A)---(10, A)---|
    // since A+C->A is very unlikely and B does not interact, we expect in the left box a reaction between
    //   7 and 8, resulting in a particle A at {0,0,-5.25}
    //   1 and 0, resulting in a particle A at {0,0,-.35}
    // in the right box we expect a reaction between
    //   4 and 5, resulting in a particle A at {0,0,1.65}
    //   9 and 10, resulting in a particle A at {0,0,5.25}
    {
        const auto particles = kernel->getKernelStateModel().getParticles();
        for (auto p : particles) readdy::log::console()->debug("particle {}", p);
        EXPECT_EQ(7, particles.size());
        EXPECT_TRUE(std::find_if(particles.begin(), particles.end(), [=](const particle_t &p) {
            return p.getType() == typeB && p.getPos() == vec_t(0, 0, -1.7);
        }) != particles.end()) << "Particle with type B should not interact with the others and thus stay where it is";
        EXPECT_TRUE(std::find_if(particles.begin(), particles.end(), [=](const particle_t &p) {
            return p.getType() == typeC && p.getPos() == vec_t(0, 0, .7);
        }) != particles.end()) << "The particle of type C is -very- unlikely to react, thus it should stay where it is";
        EXPECT_TRUE(std::find_if(particles.begin(), particles.end(), [=](const particle_t &p) {
            return p.getType() == typeC && p.getPos() == vec_t(0, 0, -1.6);
        }) != particles.end()) << "The particle of type C is -very- unlikely to react, thus it should stay where it is";
        EXPECT_TRUE(std::find_if(particles.begin(), particles.end(), [=](const particle_t &p) {
            return p.getType() == typeA && p.getPos() == vec_t(0, 0, -5.25);
        }) != particles.end()) << "This particle should be placed between the particles 7 and 8 (see above).";
        EXPECT_TRUE(std::find_if(particles.begin(), particles.end(), [=](const particle_t &p) {
            return p.getType() == typeA && p.getPos() == vec_t(0, 0, -.35);
        }) != particles.end()) << "This particle should be placed between the particles 1 and 0 (see above).";
        EXPECT_TRUE(std::find_if(particles.begin(), particles.end(), [=](const particle_t &p) {
            return p.getType() == typeA && p.getPos() == vec_t(0, 0, 1.65);
        }) != particles.end()) << "This particle should be placed between the particles 4 and 5 (see above).";
        EXPECT_TRUE(std::find_if(particles.begin(), particles.end(), [=](const particle_t &p) {
            return p.getType() == typeA && p.getPos() == vec_t(0, 0, 5.25);
        }) != particles.end()) << "This particle should be placed between the particles 9 and 10 (see above).";
    }
}


TEST(TestNextSubvolumes, ReactionRadii) {
    using fusion_t = readdy::model::reactions::Fusion;
    auto kernel = readdy::plugin::KernelProvider::getInstance().create("CPU");
    kernel->getKernelContext().setDiffusionConstant("A", 1.0);

    auto nextSubvolumes = readdy::util::static_unique_ptr_cast<readdy::kernel::cpu::actions::reactions::CPUNextSubvolumes>(
            kernel->createAction<readdy::model::actions::reactions::NextSubvolumes>(1)
    );
    kernel->getKernelContext().configure();
    EXPECT_EQ(0, nextSubvolumes->getMaxReactionRadius()) << "no reactions present, expect max radius to be 0";

    kernel->registerReaction<fusion_t>("annihilation", "A", "A", "A", 1.0, 6.0);
    kernel->getKernelContext().configure();
    EXPECT_EQ(6.0, nextSubvolumes->getMaxReactionRadius()) << "one reaction with radius 6.0";

    kernel->registerReaction<fusion_t>("annihilation2", "A", "A", "A", 1.0, 5.0);
    kernel->getKernelContext().configure();
    EXPECT_EQ(6.0, nextSubvolumes->getMaxReactionRadius()) << "max(5.0, 6.0) = 6.0";

    kernel->registerReaction<fusion_t>("annihilation3", "A", "A", "A", 1.0, 7.0);
    kernel->getKernelContext().configure();
    EXPECT_EQ(7.0, nextSubvolumes->getMaxReactionRadius()) << "max(5.0, 6.0, 7.0) = 7.0";
}
