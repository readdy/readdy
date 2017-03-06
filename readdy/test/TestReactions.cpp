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
    kernel->getKernelContext().registerParticleType("B", 2.0, 1.0);
    kernel->getKernelContext().registerParticleType("A", 1.0, 1.0);
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

    auto stdRand = [](double lower = 0.0, double upper = 1.0) -> double {
        return static_cast <double> (std::rand()) / (RAND_MAX / (upper - lower)) + lower;
    };

    kernel->getKernelContext().registerParticleType("A", 1.0, 1.0);
    kernel->getKernelContext().registerParticleType("B", 1.0, 1.0);
    kernel->getKernelContext().registerParticleType("AB", 0.0, 1.0);
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
        conf.configureAndRun(1, 10);
    }

}

INSTANTIATE_TEST_CASE_P(TestReactionsCore, TestReactions,
                        ::testing::ValuesIn(readdy::testing::getKernelsToTest()));
}