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
#include <readdy/testing/Utils.h>
#include <readdy/common/FloatingPoints.h>
#include <readdy/kernel/cpu/actions/reactions/CPUGillespie.h>
#include <readdy/model/reactions/Fusion.h>
#include <readdy/model/reactions/Fission.h>
#include <readdy/model/reactions/Decay.h>
#include <readdy/model/reactions/Enzymatic.h>
#include <readdy/model/reactions/Conversion.h>

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
    using data_t = readdy::kernel::cpu::data::DefaultDataContainer;
    auto kernel = std::make_unique<readdy::kernel::cpu::CPUKernel>();
    kernel->context().periodicBoundaryConditions() = {{false, false, false}};
    kernel->context().boxSize() = {{100, 100, 100}};
    kernel->context().particleTypes().add("A", .1); // type id 0
    kernel->context().particleTypes().add("B", .1); // type id 1
    kernel->context().particleTypes().add("C", .1); // type id 2

    // test conversion
    {
        conversion_t conversion("A->B", 0, 1, 1);
        particle_t p_A{0, 0, 0, 0};

        data_t data {kernel->context(), kernel->pool()};
        data.addParticles({p_A});

        data_t::EntriesUpdate newParticles{};
        std::vector<data_t::size_type> decayedEntries {};

        reac::performReaction(&data, kernel->context(), 0, 0, newParticles, decayedEntries, &conversion, nullptr);

        EXPECT_EQ(data.entry_at(0).type, conversion.getTypeTo());
        EXPECT_EQ(data.pos(0), readdy::Vec3(0,0,0));
    }

    // test fusion
    {
        data_t data {kernel->context(), kernel->pool()};

        data_t::EntriesUpdate newParticles{};
        std::vector<data_t::size_type> decayedEntries {};

        readdy::scalar eductDistance = .4;
        readdy::scalar weight1 = .3, weight2 = .7;
        fusion_t fusion ("A+B->C", 0, 1, 2, 1, eductDistance, weight1, weight2);
        particle_t p_A{1, 0, 0, 0};
        particle_t p_B{-1, 0, 0, 1};
        data.addParticles({p_A, p_B});

        reac::performReaction(&data, kernel->context(), 0, 1, newParticles, decayedEntries, &fusion, nullptr);
        EXPECT_EQ(decayedEntries.size(), 2);
        EXPECT_TRUE(decayedEntries.size() == 2 && (decayedEntries.at(0) == 1 || decayedEntries.at(1) == 1));
        data.update(std::make_pair(std::move(newParticles), std::move(decayedEntries)));
        EXPECT_EQ(data.entry_at(0).type, fusion.getTo());
        if(kernel->singlePrecision()) {
            EXPECT_FVEC3_EQ(readdy::Vec3(.4, 0, 0), data.pos(0));
        } else {
            EXPECT_FVEC3_EQ(readdy::Vec3(.4, 0, 0), data.pos(0));
        }
    }

    // fission
    {
        data_t data {kernel->context(), kernel->pool()};

        data_t::EntriesUpdate newParticles{};
        std::vector<data_t::size_type> decayedEntries {};

        readdy::scalar productDistance = .4;
        readdy::scalar weight1 = .3, weight2 = .7;
        fission_t fission ("C->A+B", 2, 0, 1, 1, productDistance, weight1, weight2);
        particle_t p_C{0, 0, 0, 2};
        data.addParticle(p_C);

        reac::performReaction(&data, kernel->context(), 0, 0, newParticles, decayedEntries, &fission, nullptr);
        data.update(std::make_pair(std::move(newParticles), std::move(decayedEntries)));

        EXPECT_EQ(data.entry_at(0).type, fission.getTo1());
        EXPECT_EQ(data.entry_at(1).type, fission.getTo2());
        auto p_12 = readdy::bcs::shortestDifference(data.pos(0), data.pos(1), kernel->context().boxSize(),
                                                    kernel->context().periodicBoundaryConditions());
        auto p_12_nondirect = data.pos(1) - data.pos(0);
        EXPECT_EQ(p_12_nondirect, p_12);
        auto distance = std::sqrt(p_12 * p_12);
        if(readdy::single_precision) {
            EXPECT_FLOAT_EQ(productDistance, distance);
        } else {
            EXPECT_DOUBLE_EQ(productDistance, distance);
        }
    }

    // enzymatic 1
    {
        data_t data{kernel->context(), kernel->pool()};

        data_t::EntriesUpdate newParticles{};
        std::vector<data_t::size_type > decayedEntries{};

        enzymatic_t enzymatic ("A+C->B+C", 2, 0, 1, 1, .5);
        particle_t p_A{0, 0, 0, 0};
        particle_t p_C{5, 5, 5, 2};
        data.addParticles({p_A, p_C});
        reac::performReaction(&data, kernel->context(), 0, 1, newParticles, decayedEntries, &enzymatic, nullptr);
        data.update(std::make_pair(std::move(newParticles), std::move(decayedEntries)));
        {
            const auto &e1 = data.entry_at(0);
            const auto &e2 = data.entry_at(1);
            if (e1.type == enzymatic.getCatalyst()) {
                EXPECT_EQ(enzymatic.getCatalyst(), e1.type);
                EXPECT_EQ(enzymatic.getTo(), e2.type);
                EXPECT_EQ(p_C.getPos(), e1.pos);
                EXPECT_EQ(p_A.getPos(), e2.pos);
            } else {
                EXPECT_EQ(enzymatic.getCatalyst(), e2.type);
                EXPECT_EQ(enzymatic.getTo(), e1.type);
                EXPECT_EQ(p_C.getPos(), e2.pos);
                EXPECT_EQ(p_A.getPos(), e1.pos);
            }
        }
    }
    // enzymatic 2
    {
        data_t data{kernel->context(), kernel->pool()};

        data_t::EntriesUpdate newParticles{};
        std::vector<data_t::size_type> decayedEntries{};

        enzymatic_t enzymatic ("A+C->B+C", 2, 0, 1, 1, .5);
        particle_t p_A{0, 0, 0, 0};
        particle_t p_C{5, 5, 5, 2};
        data.addParticles({p_C, p_A});
        reac::performReaction(&data, kernel->context(), 0, 1, newParticles, decayedEntries, &enzymatic, nullptr);
        data.update(std::make_pair(std::move(newParticles), std::move(decayedEntries)));
        {
            const auto &e1 = data.entry_at(0);
            const auto &e2 = data.entry_at(1);
            if (e1.type == enzymatic.getCatalyst()) {
                EXPECT_EQ(enzymatic.getCatalyst(), e1.type);
                EXPECT_EQ(enzymatic.getTo(), e2.type);
                EXPECT_EQ(p_C.getPos(), e1.pos);
                EXPECT_EQ(p_A.getPos(), e2.pos);
            } else {
                EXPECT_EQ(enzymatic.getCatalyst(), e2.type);
                EXPECT_EQ(enzymatic.getTo(), e1.type);
                EXPECT_EQ(p_C.getPos(), e2.pos);
                EXPECT_EQ(p_A.getPos(), e1.pos);
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
    kernel->context().boxSize() = {{10, 10, 10}};
    kernel->context().particleTypes().add("X", .25);
    kernel->context().reactions().addDecay("X decay", "X", 1e16);
    kernel->context().reactions().addFission("X fission", "X", "X", "X", .5, .3);

    auto &&integrator = kernel->actions().eulerBDIntegrator(1);
    auto &&forces = kernel->actions().calculateForces();
    auto &&neighborList = kernel->actions().updateNeighborList();
    auto &&reactions = kernel->actions().gillespie(1);

    auto pp_obs = kernel->observe().positions(1);
    pp_obs->callback() = ([](const readdy::model::observables::Positions::result_type &t) {
        /* ignore */
    });
    auto connection = kernel->connectObservable(pp_obs.get());

    const int n_particles = 200;
    const auto typeId = kernel->context().particleTypes().idOf("X");
    std::vector<readdy::model::Particle> particlesToBeginWith{n_particles, {0, 0, 0, typeId}};
    kernel->stateModel().addParticles(particlesToBeginWith);
    neighborList->perform();
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

TEST(CPUTestReactions, TestGillespieParallel) {

    using namespace readdy;

    using particle_t = readdy::model::Particle;
    using vec_t = readdy::Vec3;
    using fusion_t = readdy::model::reactions::Fusion;
    using fission_t = readdy::model::reactions::Fission;
    using enzymatic_t = readdy::model::reactions::Enzymatic;
    using conversion_t = readdy::model::reactions::Conversion;
    using death_t = readdy::model::reactions::Decay;
    using particle_t = readdy::model::Particle;
    auto kernel = std::make_unique<readdy::kernel::cpu::CPUKernel>();
    kernel->context().boxSize() = {{10, 10, 30}};
    kernel->context().periodicBoundaryConditions() = {{true, true, false}};

    kernel->context().particleTypes().add("A", .25);
    kernel->context().particleTypes().add("B", .25);
    kernel->context().particleTypes().add("C", .25);
    kernel->context().potentials().addBox("A", .0001, {-4.9, -4.9, -14.9}, {9.8, 9.8, 29.8});
    kernel->context().potentials().addBox("B", .0001, {-4.9, -4.9, -14.9}, {9.8, 9.8, 29.8});
    kernel->context().potentials().addBox("C", .0001, {-4.9, -4.9, -14.9}, {9.8, 9.8, 29.8});
    readdy::scalar reactionRadius = 1.0;
    kernel->context().reactions().addFusion("annihilation", "A", "A", "A", 1e16, reactionRadius);
    kernel->context().reactions().addFusion("very unlikely", "A", "C", "A", std::numeric_limits<readdy::scalar>::min(), reactionRadius);
    kernel->context().reactions().addFusion("dummy reaction", "A", "B", "A", 0.0, reactionRadius);

    const auto typeA = kernel->context().particleTypes().idOf("A");
    const auto typeB = kernel->context().particleTypes().idOf("B");
    const auto typeC = kernel->context().particleTypes().idOf("C");

    // this particle goes right into the middle, i.e., into the halo region
    kernel->stateModel().addParticle({0, 0, 0, typeA});            // 0
    // these particles go left and right of this particle into the boxes as problematic ones
    kernel->stateModel().addParticle({0, 0, static_cast<scalar>(-.7), typeA});          // 1
    kernel->stateModel().addParticle({0, 0, static_cast<scalar>(.7), typeC});           // 2
    // these particles are far enough away from the halo region but still conflict (transitively) with the 1st layer
    kernel->stateModel().addParticle({0, 0, static_cast<scalar>(-1.6), typeC});         // 3
    kernel->stateModel().addParticle({0, 0, static_cast<scalar>(1.6), typeA});          // 4
    kernel->stateModel().addParticle({0, 0, static_cast<scalar>(1.7), typeA});          // 5
    // this particle are conflicting but should not appear as their reaction rate is 0
    kernel->stateModel().addParticle({0, 0, static_cast<scalar>(-1.7), typeB});         // 6
    // these particles are well inside the boxes and should not be considered problematic
    kernel->stateModel().addParticle({0, 0, static_cast<scalar >(-5), typeA});           // 7
    kernel->stateModel().addParticle({0, 0, static_cast<scalar>(-5.5), typeA});         // 8
    kernel->stateModel().addParticle({0, 0, static_cast<scalar>(5), typeA});            // 9
    kernel->stateModel().addParticle({0, 0, static_cast<scalar>(5.5), typeA});          // 10

    // a box width in z direction of 12 should divide into two boxes of 5x5x6
    {
        fix_n_threads n_threads{kernel.get(), 2};
        EXPECT_EQ(2, kernel->getNThreads());
        auto &&neighborList = kernel->actions().updateNeighborList();
        std::unique_ptr<readdy::kernel::cpu::actions::reactions::CPUGillespie> reactions = readdy::util::static_unique_ptr_cast_no_del<readdy::kernel::cpu::actions::reactions::CPUGillespie>(
                kernel->actions().gillespie(1)
        );
        neighborList->perform();
        reactions->perform({});
        /*reactions->perform();
        EXPECT_EQ(1.0, reactions->getMaxReactionRadius());
        EXPECT_EQ(15.0, reactions->getBoxWidth());
        EXPECT_EQ(2, reactions->getLongestAxis());
        EXPECT_TRUE(reactions->getOtherAxis1() == 0 || reactions->getOtherAxis1() == 1);
        if (reactions->getOtherAxis1() == 0) {
            EXPECT_EQ(1, reactions->getOtherAxis2());
        } else {
            EXPECT_EQ(0, reactions->getOtherAxis2());
        }*/
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
        const auto particles = kernel->stateModel().getParticles();
        for (auto p : particles) readdy::log::debug("particle {}", p);
        EXPECT_EQ(7, particles.size());
        EXPECT_TRUE(std::find_if(particles.begin(), particles.end(), [=](const particle_t &p) -> bool {
            return p.getType() == typeB && p.getPos() == vec_t(0, 0, -1.7);
        }) != particles.end()) << "Particle with type B should not interact with the others and thus stay where it is";
        EXPECT_TRUE(std::find_if(particles.begin(), particles.end(), [=](const particle_t &p) -> bool {
            return p.getType() == typeC && p.getPos() == vec_t(0, 0, .7);
        }) != particles.end()) << "The particle of type C is -very- unlikely to react, thus it should stay where it is";
        EXPECT_TRUE(std::find_if(particles.begin(), particles.end(), [=](const particle_t &p) -> bool {
            return p.getType() == typeC && p.getPos() == vec_t(0, 0, -1.6);
        }) != particles.end()) << "The particle of type C is -very- unlikely to react, thus it should stay where it is";
        EXPECT_TRUE(std::find_if(particles.begin(), particles.end(), [=](const particle_t &p) -> bool {
            return p.getType() == typeA && p.getPos() == vec_t(0, 0, -5.25);
        }) != particles.end()) << "This particle should be placed between the particles 7 and 8 (see above).";
        EXPECT_TRUE(std::find_if(particles.begin(), particles.end(), [=](const particle_t &p) -> bool {
            return p.getType() == typeA && p.getPos() == vec_t(0, 0, -.35);
        }) != particles.end()) << "This particle should be placed between the particles 1 and 0 (see above).";
        auto it = std::find_if(particles.begin(), particles.end(), [=](const particle_t &p) -> bool {
            auto fpa = readdy::fp::FloatingPoint<readdy::scalar>(1.65);
            auto fpb = readdy::fp::FloatingPoint<readdy::scalar>(p.getPos().z);
            return p.getType() == typeA && fpa.AlmostEquals(fpb);
        });
        EXPECT_TRUE(it != particles.end()) << "This particle should be placed between the particles 4 and 5 (see above).";
        EXPECT_TRUE(std::find_if(particles.begin(), particles.end(), [=](const particle_t &p) -> bool {
            return p.getType() == typeA && p.getPos() == vec_t(0, 0, 5.25);
        }) != particles.end()) << "This particle should be placed between the particles 9 and 10 (see above).";
    }
}
