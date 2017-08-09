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
 * @date 21.06.16
 */

#include <readdy/plugin/KernelProvider.h>
#include <readdy/api/SimulationScheme.h>
#include <readdy/testing/KernelTest.h>
#include <readdy/testing/Utils.h>

namespace {

struct TestReactions : KernelTest {

};

TEST_P(TestReactions, TestReactionFactory) {
    kernel->getKernelContext().particle_types().add("B", 2.0, 1.0);
    kernel->getKernelContext().particle_types().add("A", 1.0, 1.0);
    kernel->registerReaction<readdy::model::reactions::Conversion>("A to B", "A", "B", 0.55);

    {
        // sanity check of operator<< for reactions
        const auto r = kernel->getReactionFactory().createReaction<readdy::model::reactions::Decay>("decay", 0, .1);
        readdy::log::debug("decay reaction: {}", *r);
    }
}

TEST_P(TestReactions, TestConstantNumberOfParticleType) {
    // scenario: two particle types A and B, which can form a complex AB which after a time is going to dissolve back
    // into A and B. Therefore, the numbers #(A) + #(AB) and #(B) + #(AB) must remain constant.

    using n_particles_obs = readdy::model::observables::NParticles;

    auto stdRand = [](readdy::scalar lower = 0.0, readdy::scalar upper = 1.0) -> readdy::scalar {
        return static_cast <readdy::scalar> (std::rand()) / (RAND_MAX / (upper - lower)) + lower;
    };

    kernel->getKernelContext().particle_types().add("A", 1.0, 1.0);
    kernel->getKernelContext().particle_types().add("B", 1.0, 1.0);
    kernel->getKernelContext().particle_types().add("AB", 0.0, 1.0);
    kernel->getKernelContext().setPeriodicBoundary(true, true, true);
    kernel->getKernelContext().setBoxSize(5, 5, 5);
    kernel->registerReaction<readdy::model::reactions::Fusion>("Form complex", "A", "B", "AB", .5, 1.0);
    kernel->registerReaction<readdy::model::reactions::Fission>("Dissolve", "AB", "A", "B", .5, 1.0);

    unsigned long n_A = 50;
    unsigned long n_B = n_A;

    for (unsigned long i = 0; i < n_A; ++i) {
        kernel->addParticle("A", {stdRand(-2.5, 2.5), stdRand(-2.5, 2.5), stdRand(-2.5, 2.5)});
        kernel->addParticle("B", {stdRand(-2.5, 2.5), stdRand(-2.5, 2.5), stdRand(-2.5, 2.5)});
    }


    auto obs = kernel->createObservable<n_particles_obs>(1, std::vector<std::string>({"A", "B", "AB"}));
    auto conn = kernel->connectObservable(obs.get());
    obs->setCallback([&n_A, &n_B](const n_particles_obs::result_t &result) {
        if (result.size() == 2) {
            EXPECT_EQ(n_A, result[0] + result[2])
                                << "Expected #(A)+#(AB)==" << n_A << ", but #(A)=" << result[0] << ", #(AB)="
                                << result[2];
            EXPECT_EQ(n_B, result[1] + result[2])
                                << "Expected #(B)+#(AB)==" << n_B << ", but #(B)=" << result[1] << ", #(AB)="
                                << result[2];
        }
    });

    {
        auto conf = readdy::api::SchemeConfigurator<readdy::api::ReaDDyScheme>(kernel.get(), true);
        const auto progs = kernel->getAvailableActions();
        if (std::find(progs.begin(), progs.end(), "GillespieParallel") != progs.end()) {
            conf = std::move(conf.withReactionScheduler<readdy::model::actions::reactions::GillespieParallel>());
        }
        conf.configureAndRun(10, 1);
    }

}

// Compare two vectors via their distance to the (1,1,1) plane in 3d space. This is quite arbitrary and only used to construct an ordered set.
struct Vec3ProjectedLess {
    bool operator()(const readdy::model::Vec3& lhs, const readdy::model::Vec3& rhs) const {
        return (lhs.x + lhs.y + lhs.z) < (rhs.x + rhs.y + rhs.z);
    }
};

TEST_P(TestReactions, FusionFissionWeights) {
    /* A diffuses, F is fixed
     * Also during reactions, weights are chosen such that F does not move.
     * Idea: position F particles and remember their positions (ordered set), do ONE time-step and check if current positions are still the same.
     */
    kernel->getKernelContext().particle_types().add("A", 0.5, 1.0);
    kernel->getKernelContext().particle_types().add("F", 0.0, 1.0);
    kernel->getKernelContext().setPeriodicBoundary(true, true, true);
    kernel->getKernelContext().setBoxSize(20, 20, 20);

    const readdy::scalar weightF {static_cast<readdy::scalar>(0)};
    const readdy::scalar weightA  {static_cast<readdy::scalar>(1.)};
    kernel->registerReaction<readdy::model::reactions::Fusion>("F+A->F", "F", "A", "F", 1.0, 2.0, weightF, weightA);
    kernel->registerReaction<readdy::model::reactions::Fission>("F->F+A", "F", "F", "A", 1.0, 2.0, weightF, weightA);

    std::set<readdy::model::Vec3, Vec3ProjectedLess> fPositions;
    auto n3 = readdy::model::rnd::normal3<readdy::scalar>;
    for (std::size_t i = 0; i < 15; ++i) {
        auto fPos = n3(static_cast<readdy::scalar>(0.), static_cast<readdy::scalar>(0.8));
        fPositions.emplace(fPos);
        kernel->addParticle("F", fPos);
        kernel->addParticle("A", n3(static_cast<readdy::scalar>(0.), static_cast<readdy::scalar>(1.)));
    }

    auto obs = kernel->createObservable<readdy::model::observables::Positions>(1, std::vector<std::string>({"F"}));
    obs->setCallback(
            [&fPositions](const readdy::model::observables::Positions::result_t &result) {
                std::set<readdy::model::Vec3, Vec3ProjectedLess> checklist;
                for (const auto &pos : result) {
                    checklist.emplace(pos);
                }
                EXPECT_EQ(fPositions.size(), checklist.size());
                auto itPos = fPositions.begin();
                auto itCheck = checklist.begin();
                while (itPos != fPositions.end()) {
                    EXPECT_VEC3_NEAR(*itPos, *itCheck, readdy::double_precision ? 1e-8 : 1e-6);
                    ++itPos;
                    ++itCheck;
                }
            }
    );
    auto connection = kernel->connectObservable(obs.get());

    {
        auto conf = readdy::api::SchemeConfigurator<readdy::api::ReaDDyScheme>(kernel.get(), true);
        const auto progs = kernel->getAvailableActions();
        if (std::find(progs.begin(), progs.end(), "GillespieParallel") != progs.end()) {
            conf = std::move(conf.withReactionScheduler<readdy::model::actions::reactions::GillespieParallel>());
        }
        conf.configureAndRun(1, 0.5);
    }
}

/*
 * @todo this is rather an integration test that should be separated from the rest
 * TEST_P(TestReactions, ConstantNumberOfParticles) {
 *
    using namespace readdy;
    // A is absorbed and created by F, while the number of F stays constant, this test spans multiple timesteps
    kernel->getKernelContext().particle_types().add("A", 0.5, 1.0);
    kernel->getKernelContext().particle_types().add("F", 0.0, 1.0);
    kernel->getKernelContext().setPeriodicBoundary(true, true, true);
    kernel->getKernelContext().setBoxSize(10, 10, 10);

    const readdy::scalar weightF = 0.;
    const readdy::scalar weightA = 1.;
    kernel->registerReaction<readdy::model::reactions::Fusion>("F+A->F", "F", "A", "F", .1, 2.0, weightF, weightA);

    auto n3 = readdy::model::rnd::normal3<>;
    for (std::size_t i = 0; i < 200; ++i) {
        kernel->addParticle("F", n3(0., 1.));
        kernel->addParticle("A", n3(0., 1.));
    }

    auto obs = kernel->createObservable<readdy::model::observables::NParticles>(1, std::vector<std::string>({"F"}));
    obs->setCallback(
            [](const readdy::model::observables::NParticles::result_t &result) {
                EXPECT_EQ(result[0], 200);
                log::error("::: {}", result[0]);
            }
    );
    auto connection = kernel->connectObservable(obs.get());

    {
        auto conf = readdy::api::SchemeConfigurator<readdy::api::ReaDDyScheme>(kernel.get(), true);
        conf.withReactionScheduler<readdy::model::actions::reactions::GillespieParallel>();
        conf.withSkinSize(5.);
        conf.configureAndRun(400, 0.01);
    }
}*/

INSTANTIATE_TEST_CASE_P(TestReactionsCore, TestReactions,
                        ::testing::ValuesIn(readdy::testing::getKernelsToTest()));
}