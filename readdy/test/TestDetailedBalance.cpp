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
 * @file TestDetailedBalance.cpp
 * @brief Kernel-non-specific tests of the detailed-balance reaction handler
 * @author chrisfroe
 * @date 29.05.18
 * @copyright BSD-3
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <readdy/common/common.h>
#include <readdy/testing/KernelTest.h>
#include <readdy/testing/Utils.h>
#include <readdy/model/actions/Action.h>

namespace m = readdy::model;
namespace r = m::actions::reactions;

using namespace readdytesting::kernel;

TEST_CASE("Test detailed balance config.", "[detailed-balance]") {
    m::Context ctx;
    ctx.kBT() = 1;
    ctx.particleTypes().add("A", 1.);
    ctx.particleTypes().add("B", 1.);
    ctx.particleTypes().add("C", 1.);

    SECTION("Config values fusion and fission") {
        ctx.boxSize() = {21.544346900318832, 21.544346900318832, 21.544346900318832};

        ctx.potentials().addHarmonicRepulsion("A", "B", 2., 3.);
        ctx.reactions().add("fusion: A +(5) B -> C", 1.);
        ctx.reactions().add("fission: C -> A +(5) B", 2.);

        auto idFusion = ctx.reactions().idOf("fusion");
        auto idFission = ctx.reactions().idOf("fission");
        r::ReversibleReactionConfig conf(idFusion, idFission, ctx);
        REQUIRE(conf.totalVolume == Catch::Approx(10000).epsilon(1e-4));
        REQUIRE(conf.kbt == 1.);
        REQUIRE(conf.lhsInteractionRadius == 3.);
        REQUIRE(conf.reactionRadius == 5.);
        REQUIRE(conf.numberLhsTypes == 2);
        REQUIRE(conf.numberRhsTypes == 1);
        REQUIRE(conf.rhsTypes[0] == ctx.particleTypes().idOf("C"));
        REQUIRE(conf.lhsTypes[0] == ctx.particleTypes().idOf("A"));
        REQUIRE(conf.lhsTypes[1] == ctx.particleTypes().idOf("B"));
        REQUIRE(conf.lhsInteractionVolume == Catch::Approx(113.09733).epsilon(1e-4));
        REQUIRE(conf.effectiveLhsInteractionVolume == Catch::Approx(68.09910).epsilon(1e-4));
        REQUIRE(conf.effectiveLhsReactionVolume == Catch::Approx(478.60057).epsilon(1e-4));
        REQUIRE(conf.equilibriumConstant == Catch::Approx(41.60045).epsilon(1e-4));
        REQUIRE(conf.macroBackwardRate == Catch::Approx(2.));
        REQUIRE(conf.macroForwardRate == Catch::Approx(480.76392).epsilon(1e-4));
        REQUIRE(conf.acceptancePrefactor == 1.);
    }

    SECTION("Config values conversion and conversion") {
        ctx.boxSize() = {10., 10., 10.};
        ctx.reactions().add("forward: A -> B", 1.);
        ctx.reactions().add("backward: B -> A", 2.);

        auto idForward = ctx.reactions().idOf("forward");
        auto idBackward = ctx.reactions().idOf("backward");

        r::ReversibleReactionConfig conf(idForward, idBackward, ctx);
        REQUIRE(conf.totalVolume == Catch::Approx(1000).epsilon(1e-4));
        REQUIRE(conf.numberLhsTypes == 1);
        REQUIRE(conf.numberRhsTypes == 1);
        REQUIRE(conf.lhsTypes[0] == ctx.particleTypes().idOf("A"));
        REQUIRE(conf.rhsTypes[0] == ctx.particleTypes().idOf("B"));
        REQUIRE(conf.lhsInteractionRadius == 0.);
        REQUIRE(conf.lhsInteractionVolume == 0.);
        REQUIRE(conf.rhsInteractionRadius == 0.);
        REQUIRE(conf.rhsInteractionVolume == 0.);

        REQUIRE(conf.effectiveLhsInteractionVolume == 0.);
        REQUIRE(conf.effectiveLhsReactionVolume == 0.);
        REQUIRE(conf.effectiveRhsInteractionVolume == 0.);
        REQUIRE(conf.effectiveRhsReactionVolume == 0.);

        REQUIRE(conf.acceptancePrefactor == 1.);

        REQUIRE(conf.macroForwardRate == 1.);
        REQUIRE(conf.macroBackwardRate == 2.);
        REQUIRE(conf.equilibriumConstant == 2.);
    }

    SECTION("Config values enzymatic and enzymatic") {
        ctx.boxSize() = {10., 10., 10.};
        ctx.potentials().addHarmonicRepulsion("A", "C", 2., 3.);
        ctx.potentials().addHarmonicRepulsion("B", "C", 2., 3.);
        ctx.reactions().add("forward: A +(5) C -> B + C", 1.);
        ctx.reactions().add("backward: B +(5) C -> A + C", 2.);

        auto idForward = ctx.reactions().idOf("forward");
        auto idBackward = ctx.reactions().idOf("backward");
        r::ReversibleReactionConfig conf(idForward, idBackward, ctx);
        REQUIRE(conf.totalVolume == Catch::Approx(1000).epsilon(1e-4));
        REQUIRE(conf.numberLhsTypes == 2);
        REQUIRE(conf.numberRhsTypes == 2);
        REQUIRE(conf.lhsTypes[0] == ctx.particleTypes().idOf("A"));
        REQUIRE(conf.rhsTypes[0] == ctx.particleTypes().idOf("B"));
        REQUIRE(conf.lhsTypes[1] == ctx.particleTypes().idOf("C"));
        REQUIRE(conf.rhsTypes[1] == ctx.particleTypes().idOf("C"));
        REQUIRE(conf.lhsInteractionRadius == 3.);
        REQUIRE(conf.rhsInteractionRadius == 3.);
        auto expectedVolume = 4. / 3. * readdy::util::numeric::pi<readdy::scalar>() * std::pow(3., 3);
        REQUIRE(conf.lhsInteractionVolume == Catch::Approx(expectedVolume));
        REQUIRE(conf.rhsInteractionVolume == Catch::Approx(expectedVolume));

        REQUIRE(conf.effectiveLhsInteractionVolume == Catch::Approx(68.09910).epsilon(1e-4));
        REQUIRE(conf.effectiveLhsReactionVolume == Catch::Approx(478.60057).epsilon(1e-4));
        REQUIRE(conf.effectiveRhsInteractionVolume == Catch::Approx(68.09910).epsilon(1e-4));
        REQUIRE(conf.effectiveRhsReactionVolume == Catch::Approx(478.60057).epsilon(1e-4));

        INFO(fmt::format("ratio of effective reaction volumina is expected to be 1., got {}", conf.acceptancePrefactor));
        REQUIRE(conf.acceptancePrefactor == Catch::Approx(1.).epsilon(1e-4));

        REQUIRE(conf.macroForwardRate == Catch::Approx(501.15150).epsilon(1e-4));
        REQUIRE(conf.macroBackwardRate == Catch::Approx(1002.30300).epsilon(1e-4));
        REQUIRE(conf.equilibriumConstant == 2.);
    }

    SECTION("Invalid reaction configuration") {
        SECTION("Radii do not match"){
            ctx.reactions().add("forward: A +(4) C -> B + C", 1.);
            ctx.reactions().add("backward: B +(5) C -> A + C", 2.);
            auto idForward = ctx.reactions().idOf("forward");
            auto idBackward = ctx.reactions().idOf("backward");
            REQUIRE_THROWS_AS(r::ReversibleReactionConfig (idForward, idBackward, ctx), std::logic_error);
        }

        SECTION("Not reversible wrt species Conversion+Conversion"){
            ctx.reactions().add("forward: A -> B", 1.);
            ctx.reactions().add("backward: B -> C", 2.);
            auto idForward = ctx.reactions().idOf("forward");
            auto idBackward = ctx.reactions().idOf("backward");
            REQUIRE_THROWS_AS(r::ReversibleReactionConfig (idForward, idBackward, ctx), std::logic_error);
        }

        SECTION("Not reversible wrt species Fusion+Fission"){
            ctx.reactions().add("forward: A +(2) B -> C", 1.);
            ctx.reactions().add("backward: C -> A +(2) A", 2.);
            auto idForward = ctx.reactions().idOf("forward");
            auto idBackward = ctx.reactions().idOf("backward");
            REQUIRE_THROWS_AS(r::ReversibleReactionConfig(idForward, idBackward, ctx), std::logic_error);
        }

        SECTION("Not reversible wrt reaction radii"){
            ctx.reactions().add("forward: A +(2) B -> C", 1.);
            ctx.reactions().add("backward: C -> A +(3) B", 2.);
            auto idForward = ctx.reactions().idOf("forward");
            auto idBackward = ctx.reactions().idOf("backward");
            REQUIRE_THROWS_AS(r::ReversibleReactionConfig(idForward, idBackward, ctx), std::logic_error);
        }
    }

    SECTION("Valid reaction configuration") {
        SECTION("Fusion+Fission can be swapped") {
            ctx.reactions().add("forward: C -> B +(2) A", 2.);
            ctx.reactions().add("backward: A +(2) B -> C", 1.);
            auto idForward = ctx.reactions().idOf("forward");
            auto idBackward = ctx.reactions().idOf("backward");
            REQUIRE_NOTHROW(r::ReversibleReactionConfig (idForward, idBackward, ctx));
        }
        SECTION("Enzymatic+Enzymatic without potentials") {
            // EnzymaticEnzymatic without potentials -> acceptancePrefactor = 1
            ctx.reactions().add("forward: A +(4) C -> B + C", 1.);
            ctx.reactions().add("backward: B +(4) C -> A + C", 2.);
            auto idForward = ctx.reactions().idOf("forward");
            auto idBackward = ctx.reactions().idOf("backward");
            r::ReversibleReactionConfig conf(idForward, idBackward, ctx);
            REQUIRE(conf.acceptancePrefactor == 1.);
        }
    }
}

TEMPLATE_TEST_CASE("Test detailed balance action.", "[detailed-balance]", SingleCPU) {
    auto kernel = create<TestType>();
    auto &ctx = kernel->context();
    ctx.kBT() = 1;
    ctx.particleTypes().add("A", 1.);
    ctx.particleTypes().add("B", 1.);
    ctx.particleTypes().add("C", 1.);

    SECTION("Search reversible reactions") {
        SECTION("Enzymatic+Enzymatic") {
            ctx.reactions().add("forward: A +(4) C -> B + C", 1.);
            ctx.reactions().add("backward: B +(4) C -> A + C", 2.);
            auto &&reactions = kernel->actions().detailedBalance(1.);
            const auto &reversibleReactions = reactions->reversibleReactions();
            REQUIRE(reversibleReactions.size() == 1);
            auto reversibleType = reversibleReactions.front()->reversibleType;
            REQUIRE(reversibleType == r::ReversibleType::EnzymaticEnzymatic);
        }
        SECTION("Fusion+Fission") {
            ctx.reactions().add("fusion: A +(2) B -> C", 1.);
            ctx.reactions().add("fission: C -> A +(2) B", 2.);
            auto &&reactions = kernel->actions().detailedBalance(1.);
            const auto &reversibleReactions = reactions->reversibleReactions();
            REQUIRE(reversibleReactions.size() == 1);
            REQUIRE(reversibleReactions.front()->reversibleType == r::ReversibleType::FusionFission);
        }
        SECTION("Fission+Fusion") {
            ctx.reactions().add("fission: C -> A +(2) B", 2.);
            ctx.reactions().add("fusion: B +(2) A -> C", 1.);
            auto &&reactions = kernel->actions().detailedBalance(1.);
            const auto &reversibleReactions = reactions->reversibleReactions();
            REQUIRE(reversibleReactions.size() == 1);
            REQUIRE(reversibleReactions.front()->reversibleType == r::ReversibleType::FusionFission);
        }
        SECTION("Michaelis-Menten") {
            // one should be able to have a chain of reversible reactions like in Michaelis Menten type
            // reaction A + B <--> C <--> A + D
            // add a non-reversible, B --> D
            ctx.particleTypes().add("D", 1.);
            ctx.reactions().add("fusion1: A +(2) B -> C", 2.);
            ctx.reactions().add("fission1: C -> A +(2) B", 1.);
            ctx.reactions().add("fission2: C -> A +(3) D", 3.);
            ctx.reactions().add("fusion2: A +(3) D -> C", 4.);
            ctx.reactions().add("shortcut: B -> D", 5.);

            auto &&reactions = kernel->actions().detailedBalance(1.);
            const auto &reversibleReactions = reactions->reversibleReactions();
            REQUIRE(reversibleReactions.size() == 2);
            REQUIRE(reversibleReactions[0]->reversibleType == r::ReversibleType::FusionFission);
            REQUIRE(reversibleReactions[1]->reversibleType == r::ReversibleType::FusionFission);
            std::size_t index1, index2;
            if (reversibleReactions[0]->forwardName == "fusion1") {
                index1 = 0;
                index2 = 1;
            } else {
                index1 = 1;
                index2 = 0;
            }
            REQUIRE(reversibleReactions[index1]->forwardName == "fusion1");
            REQUIRE(reversibleReactions[index2]->forwardName == "fusion2");
            REQUIRE(reversibleReactions[index1]->backwardName == "fission1");
            REQUIRE(reversibleReactions[index2]->backwardName == "fission2");
        }
        SECTION("Many mixed") {
            // A + B <--> C, reversible
            // A + C --> B, non reversible
            // A + B <--> C + B, reversible, B is catalyst
            // B + A --> C + A, non-reversible, A is catalyst
            // A <--> B, reversible
            // C --> B, non-reversible
            ctx.reactions().add("fusionRev: A +(2) B -> C", 2.);
            ctx.reactions().add("fissionRev: C -> A +(2) B", 1.);
            ctx.reactions().add("fusionNonRev: A +(3) C -> B", 4.);
            ctx.reactions().add("enzRev1: A +(2) B -> C + B", 3.);
            ctx.reactions().add("enzRev2: C +(2) B -> A + B", 5.);
            ctx.reactions().add("enzNonRev: B +(4) A -> C + A", 6.);
            ctx.reactions().add("convRev1: A -> B", 7.);
            ctx.reactions().add("convRev2: B -> A", 8.);
            ctx.reactions().add("convNonRev: C -> B", 9.);

            auto &&reactions = kernel->actions().detailedBalance(1.);
            const auto &reversibleReactions = reactions->reversibleReactions();
            REQUIRE(reversibleReactions.size() == 3);
            for (const auto &rev : reversibleReactions) {
                switch (rev->reversibleType) {
                    case readdy::model::actions::reactions::FusionFission: {
                        REQUIRE(rev->forwardName == "fusionRev");
                        REQUIRE(rev->backwardName == "fissionRev");
                        break;
                    }
                    case readdy::model::actions::reactions::ConversionConversion: {
                        bool convrev12 = (rev->forwardName == "convRev1" and rev->backwardName == "convRev2");
                        bool convrev21 = (rev->forwardName == "convRev2" and rev->backwardName == "convRev1");
                        REQUIRE((convrev12 || convrev21));
                        break;
                    }
                    case readdy::model::actions::reactions::EnzymaticEnzymatic: {
                        bool enzrev12 = (rev->forwardName == "enzRev1" and rev->backwardName == "enzRev2");
                        bool enzrev21 = (rev->forwardName == "enzRev2" and rev->backwardName == "enzRev1");
                        REQUIRE((enzrev12 || enzrev21));
                        break;
                    }
                    default:
                        FAIL("Each type should be here");
                }
            }
        }
    }
    SECTION("Duplicate reactions") {
        ctx.reactions().add("fusion: A +(2) B -> C", 2.);
        ctx.reactions().add("fission: C -> A +(2) B", 1.);
        ctx.reactions().add("anotherFusion: A +(2) B -> C", 2.);

        REQUIRE_THROWS_AS(kernel->actions().detailedBalance(1.), std::logic_error);
    }
    SECTION("Wrong weights") {
        ctx.reactions().addFusion("fusion", "A", "B", "C", 1., 2., 0.5, 0.5);
        ctx.reactions().addFission("fission", "C", "A", "B", 1., 2., 0.3, 0.7);

        auto &&reactions = kernel->actions().detailedBalance(1.);
        const auto &reversibleReactions = reactions->reversibleReactions();
        REQUIRE(reversibleReactions.empty());
    }
}

void perform(readdy::model::Kernel *kernel, size_t nSteps, readdy::scalar timeStep, bool withIntegrator = false) {
    auto &&integrator = kernel->actions().eulerBDIntegrator(timeStep);
    auto &&forces = kernel->actions().calculateForces();
    using update_nl = readdy::model::actions::CreateNeighborList;
    auto &&initNeighborList = kernel->actions().createNeighborList(kernel->context().calculateMaxCutoff());
    auto &&neighborList = kernel->actions().updateNeighborList();
    auto &&reactions = kernel->actions().detailedBalance(timeStep);

    initNeighborList->perform();
    neighborList->perform();
    forces->perform();
    kernel->evaluateObservables(0);
    for (size_t t = 1; t < nSteps + 1; t++) {
        if (withIntegrator) {
            integrator->perform();
        }
        neighborList->perform();
        forces->perform();
        reactions->perform();
        neighborList->perform();
        forces->perform();
        kernel->evaluateObservables(t);
    }
}

TEMPLATE_TEST_CASE("Detailed balance integration tests.", "[detailed-balance]", SingleCPU) {
    auto kernel = create<TestType>();
    auto &ctx = kernel->context();
    ctx.boxSize() = {{12, 12, 12}};
    ctx.particleTypes().add("A", 1.);
    ctx.particleTypes().add("B", 1.);
    ctx.particleTypes().add("C", 1.);
    readdy::scalar reactionRadius = 2.;
    ctx.potentials().addHarmonicRepulsion("A", "B", 10., reactionRadius);
    ctx.potentials().addHarmonicRepulsion("B", "B", 10., reactionRadius);
    ctx.potentials().addHarmonicRepulsion("A", "A", 10., reactionRadius);
    ctx.potentials().addHarmonicRepulsion("A", "C", 10., 1 + 1.2599210498948732);
    ctx.potentials().addHarmonicRepulsion("B", "C", 10., 1 + 1.2599210498948732);
    ctx.potentials().addHarmonicRepulsion("C", "C", 10., 2. * 1.2599210498948732);

    SECTION("Conversion") {
        SECTION("Should be rejected") {
            ctx.reactions().addConversion("convAB", "A", "B", 1e15);
            ctx.reactions().addConversion("convBA", "B", "A", 1e-15);
            ctx.particleTypes().add("D", 0.);
            ctx.potentials().addHarmonicRepulsion("B", "D", 10., 2.);

            const auto idForward = ctx.reactions().idOf("convAB");
            const auto idBackward = ctx.reactions().idOf("convBA");

            auto countsObs = kernel->observe().reactionCounts(1, [&idForward, &idBackward](const readdy::model::observables::ReactionCounts::result_type &result) {
                const auto& counts = std::get<0>(result);
                if (counts.empty()) {
                    readdy::log::trace("reaction counts is empty, no reaction handler ran so far, skip test");
                    return;
                }
                REQUIRE(counts.at(idBackward) == 0);
                REQUIRE(counts.at(idForward) == 0);
            });
            auto countsConnection = kernel->connectObservable(countsObs.get());

            std::vector<std::string> typesToCount = {"A", "B", "D"};
            auto numbersObs = kernel->observe().nParticles(1, typesToCount, [](const readdy::model::observables::NParticles::result_type &result) {
                REQUIRE(result[0] + result[1] == 1); // conservation of A + B
                REQUIRE(result[2] == 500); // conservation of D
            });
            auto numbersConnection = kernel->connectObservable(numbersObs.get());

            const auto ida = ctx.particleTypes().idOf("A");
            const auto idd = ctx.particleTypes().idOf("D");
            kernel->stateModel().addParticle({{0., 0., 0.}, ida});

            // induce rejection by many repulsing particles in the vicinity
            for (std::size_t i = 0; i < 500; ++i) {
                auto pos = m::rnd::normal3(0., 1.);
                kernel->stateModel().addParticle({pos, idd});
            }

            readdy::scalar timeStep = 0.1;
            perform(kernel.get(), 1, timeStep);
        }

        SECTION("Should be accepted") {
            ctx.reactions().addConversion("convAB", "A", "B", 1e15);
            ctx.reactions().addConversion("convBA", "B", "A", 1e-15);
            ctx.particleTypes().add("D", 0.);
            ctx.potentials().addHarmonicRepulsion("A", "D", 10., 2.); // A and D interact, thus the initial state is unfavorable

            const auto idForward = ctx.reactions().idOf("convAB");
            const auto idBackward = ctx.reactions().idOf("convBA");

            auto countsObs = kernel->observe().reactionCounts(1,[&idForward, &idBackward](const readdy::model::observables::ReactionCounts::result_type &result) {
                const auto& counts = std::get<0>(result);
                if (counts.empty()) {
                    readdy::log::trace("reaction counts is empty, no reaction handler ran so far, skip test");
                    return;
                }
                REQUIRE(counts.at(idBackward) == 0);
                REQUIRE(counts.at(idForward) == 1);
            });
            auto countsConnection = kernel->connectObservable(countsObs.get());

            std::vector<std::string> typesToCount = {"A", "B", "D"};
            auto numbersObs = kernel->observe().nParticles(1, typesToCount, [](const readdy::model::observables::NParticles::result_type &result) {
                REQUIRE(result[0] + result[1] == 1); // conservation of A + B
                REQUIRE(result[2] == 500); // conservation of D
            });
            auto numbersConnection = kernel->connectObservable(numbersObs.get());

            const auto ida = ctx.particleTypes().idOf("A");
            const auto idd = ctx.particleTypes().idOf("D");
            kernel->stateModel().addParticle({{0., 0., 0.}, ida});

            // induce rejection by many repulsing particles in the vicinity
            for (std::size_t i = 0; i < 500; ++i) {
                auto pos = m::rnd::normal3(0., 1.);
                kernel->stateModel().addParticle({pos, idd});
            }

            readdy::scalar timeStep = 0.1;
            perform(kernel.get(), 1, timeStep);
        }
    }

    SECTION("Fusion") {
        SECTION("Should be accepted") {
            // high on rate, small off rate
            ctx.reactions().addFusion("fusion", "A", "B", "C", 1e15, reactionRadius);
            ctx.reactions().addFission("fission", "C", "A", "B", 1e-15, reactionRadius);
            ctx.particleTypes().add("D", 0.);
            ctx.potentials().addHarmonicRepulsion("A", "D", 10., 2.);
            ctx.potentials().addHarmonicRepulsion("B", "D", 10., 2.);

            const auto idfus = ctx.reactions().idOf("fusion");
            const auto idfis = ctx.reactions().idOf("fission");

            auto countsObs = kernel->observe().reactionCounts(1,[&idfus, &idfis](const readdy::model::observables::ReactionCounts::result_type &result) {
                const auto& counts = std::get<0>(result);
                if (counts.empty()) {
                    readdy::log::trace("reaction counts is empty, no reaction handler ran so far, skip test");
                    return;
                }
                REQUIRE(counts.at(idfus) == 1); // fusion must occur
                REQUIRE(counts.at(idfis) == 0); // fission shall not occur
            });
            auto countsConnection = kernel->connectObservable(countsObs.get());

            std::vector<std::string> typesToCount = {"A", "B", "C", "D"};
            auto numbersObs = kernel->observe().nParticles(1, typesToCount, [](const readdy::model::observables::NParticles::result_type &result) {
                REQUIRE(result[0] + result[2] == 1); // conservation of A + C
                REQUIRE(result[1] + result[2] == 1); // conservation of B + C
                REQUIRE(result[3] == 500); // conservation of D
                if (result[0] + result[2] != 1) {
                    readdy::log::trace("A {} B {} C {}", result[0], result[1], result[2]);
                }
            });
            auto numbersConnection = kernel->connectObservable(numbersObs.get());

            const auto ida = ctx.particleTypes().idOf("A");
            const auto idb = ctx.particleTypes().idOf("B");
            const auto idd = ctx.particleTypes().idOf("D");
            kernel->stateModel().addParticle({{0.1, 0.1, 0.1}, ida});
            kernel->stateModel().addParticle({{-0.1, -0.1, -0.1}, idb});

            // induce rejection by many repulsing particles in the vicinity interacting with the A and B
            for (std::size_t i = 0; i < 500; ++i) {
                auto pos = m::rnd::normal3(0., 1.);
                kernel->stateModel().addParticle({pos, idd});
            }

            readdy::scalar timeStep = 0.1;
            perform(kernel.get(), 1, timeStep);
        }
        SECTION("Should be rejected") {
            ctx.reactions().addFusion("fusion", "A", "B", "C", 1e15, reactionRadius);
            ctx.reactions().addFission("fission", "C", "A", "B", 1e-15, reactionRadius);
            ctx.particleTypes().add("D", 0.);
            ctx.potentials().addHarmonicRepulsion("C", "D", 10., 2.);

            const auto idfus = ctx.reactions().idOf("fusion");
            const auto idfis = ctx.reactions().idOf("fission");

            auto countsObs = kernel->observe().reactionCounts(1,[&idfus, &idfis](const readdy::model::observables::ReactionCounts::result_type &result) {
                const auto& counts = std::get<0>(result);
                if (counts.empty()) {
                    readdy::log::trace("reaction counts is empty, no reaction handler ran so far, skip test");
                    return;
                }
                // fusion shall not occur
                REQUIRE(counts.at(idfus) == 0);
                // fission shall not occur
                REQUIRE(counts.at(idfis) == 0);
            });
            auto countsConnection = kernel->connectObservable(countsObs.get());

            std::vector<std::string> typesToCount = {"A", "B", "C", "D"};
            auto numbersObs = kernel->observe().nParticles(1, typesToCount, [](const readdy::model::observables::NParticles::result_type &result) {
                // conservation of A + C
                REQUIRE(result[0] + result[2] == 1);
                // conservation of B + C
                REQUIRE(result[1] + result[2] == 1);
                // conservation of D
                REQUIRE(result[3] == 500);
                if (result[0] + result[2] != 1) {
                    readdy::log::trace("A {} B {} C {}", result[0], result[1], result[2]);
                }
            });
            auto numbersConnection = kernel->connectObservable(numbersObs.get());

            const auto ida = ctx.particleTypes().idOf("A");
            const auto idb = ctx.particleTypes().idOf("B");
            const auto idd = ctx.particleTypes().idOf("D");
            kernel->stateModel().addParticle({{0.1, 0.1, 0.1}, ida});
            kernel->stateModel().addParticle({{-0.1, -0.1, -0.1}, idb});

            // induce rejection by many repulsing particles in the vicinity
            for (std::size_t i = 0; i < 500; ++i) {
                auto pos = m::rnd::normal3(0., 1.);
                kernel->stateModel().addParticle({pos, idd});
            }

            readdy::scalar timeStep = 0.1;
            perform(kernel.get(), 1, timeStep);
        }
    }

    SECTION("Fission") {
        SECTION("Should be accepted") {
            // very high off rate
            ctx.reactions().addFusion("fusion", "A", "B", "C", 1e-15, reactionRadius);
            ctx.reactions().addFission("fission", "C", "A", "B", 1e15, reactionRadius);
            ctx.particleTypes().add("D", 0.);
            ctx.potentials().addHarmonicRepulsion("C", "D", 10., 2.);

            const auto idfus = ctx.reactions().idOf("fusion");
            const auto idfis = ctx.reactions().idOf("fission");

            auto countsObs = kernel->observe().reactionCounts(1, [&idfus, &idfis](const readdy::model::observables::ReactionCounts::result_type &result) {
                const auto& counts = std::get<0>(result);
                if (counts.empty()) {
                    readdy::log::trace("reaction counts is empty, no reaction handler ran so far, skip test");
                    return;
                }
                REQUIRE(counts.at(idfis) == 1); // fission must occur
                REQUIRE(counts.at(idfus) == 0); // fusion shall not occur
            });
            auto countsConnection = kernel->connectObservable(countsObs.get());

            std::vector<std::string> typesToCount = {"A", "B", "C", "D"};
            auto numbersObs = kernel->observe().nParticles(1, typesToCount, [](const readdy::model::observables::NParticles::result_type &result) {
                REQUIRE(result[0] + result[2] == 1); // conservation of A + C
                REQUIRE(result[1] + result[2] == 1); // conservation of B + C
                REQUIRE(result[3] == 500); // conservation of D
                if (result[0] + result[2] != 1) {
                    readdy::log::trace("A {} B {} C {}", result[0], result[1], result[2]);
                }
            });
            auto numbersConnection = kernel->connectObservable(numbersObs.get());

            const auto idc = ctx.particleTypes().idOf("C");
            const auto idd = ctx.particleTypes().idOf("D");
            kernel->stateModel().addParticle({{-0.1, -0.1, -0.1}, idc});

            // induce rejection by many repulsing particles in the vicinity
            for (std::size_t i = 0; i < 500; ++i) {
                auto pos = m::rnd::normal3(0., 1.);
                kernel->stateModel().addParticle({pos, idd});
            }

            readdy::scalar timeStep = 0.1;
            perform(kernel.get(), 1, timeStep);
        }
        SECTION("Should be rejected") {
            ctx.reactions().addFusion("fusion", "A", "B", "C", 1e-15, reactionRadius);
            ctx.reactions().addFission("fission", "C", "A", "B", 1e15, reactionRadius);
            ctx.particleTypes().add("D", 0.);
            ctx.potentials().addHarmonicRepulsion("A", "D", 10., 2.);
            ctx.potentials().addHarmonicRepulsion("B", "D", 10., 2.);

            const auto idfus = ctx.reactions().idOf("fusion");
            const auto idfis = ctx.reactions().idOf("fission");

            auto countsObs = kernel->observe().reactionCounts(1,[&idfus, &idfis](const readdy::model::observables::ReactionCounts::result_type &result) {
                const auto& counts = std::get<0>(result);
                if (counts.empty()) {
                    readdy::log::trace("reaction counts is empty, no reaction handler ran so far, skip test");
                    return;
                }
                REQUIRE(counts.at(idfis) == 0); // fission shall not occur
                REQUIRE(counts.at(idfus) == 0); // fusion shall not occur
            });
            auto countsConnection = kernel->connectObservable(countsObs.get());

            std::vector<std::string> typesToCount = {"A", "B", "C", "D"};
            auto numbersObs = kernel->observe().nParticles(1, typesToCount, [](const readdy::model::observables::NParticles::result_type &result) {
                REQUIRE(result[0] + result[2] == 1); // conservation of A + C
                REQUIRE(result[1] + result[2] == 1); // conservation of B + C
                REQUIRE(result[3] == 500); // conservation of D
                if (result[0] + result[2] != 1) {
                    readdy::log::trace("A {} B {} C {}", result[0], result[1], result[2]);
                }
            });
            auto numbersConnection = kernel->connectObservable(numbersObs.get());

            const auto idc = ctx.particleTypes().idOf("C");
            const auto idd = ctx.particleTypes().idOf("D");
            kernel->stateModel().addParticle({{-0.1, -0.1, -0.1}, idc});

            // induce rejection by many repulsing particles in the vicinity
            for (std::size_t i = 0; i < 500; ++i) {
                auto pos = m::rnd::normal3(0., 1.);
                kernel->stateModel().addParticle({pos, idd});
            }

            readdy::scalar timeStep = 0.1;
            perform(kernel.get(), 1, timeStep);
        }
    }

    SECTION("Enzymatic") {
        readdy::scalar interactionDistance = 2.;
        SECTION("Should be accepted") {
            ctx.reactions().addEnzymatic("convAB", "C", "A", "B", 1e15, interactionDistance);
            ctx.reactions().addEnzymatic("convBA", "C", "B", "A", 1e-15, interactionDistance);
            ctx.particleTypes().add("D", 0.);
            ctx.potentials().addHarmonicRepulsion("A", "D", 10., 2.); // A and D interact, thus the initial state is unfavorable

            const auto idForward = ctx.reactions().idOf("convAB");
            const auto idBackward = ctx.reactions().idOf("convBA");

            auto countsObs = kernel->observe().reactionCounts(1, [&idForward, &idBackward](const readdy::model::observables::ReactionCounts::result_type &result) {
                const auto& counts = std::get<0>(result);
                if (counts.empty()) {
                    readdy::log::trace("reaction counts is empty, no reaction handler ran so far, skip test");
                    return;
                }
                REQUIRE(counts.at(idBackward) == 0);
                REQUIRE(counts.at(idForward) == 1);
            });
            auto countsConnection = kernel->connectObservable(countsObs.get());

            std::vector<std::string> typesToCount = {"A", "B", "C", "D"};
            auto numbersObs = kernel->observe().nParticles(1, typesToCount, [](const readdy::model::observables::NParticles::result_type &result) {
                REQUIRE(result[0] + result[1] == 1); // conservation of A + B
                REQUIRE(result[2] == 1); // conservation of C
                REQUIRE(result[3] == 500); // conservation of D
            });
            auto numbersConnection = kernel->connectObservable(numbersObs.get());

            const auto ida = ctx.particleTypes().idOf("A");
            const auto idc = ctx.particleTypes().idOf("C");
            const auto idd = ctx.particleTypes().idOf("D");
            kernel->stateModel().addParticle({{0., 0., 0.}, ida});
            kernel->stateModel().addParticle({{0., 0., 1.}, idc});

            // induce rejection by many repulsing particles in the vicinity
            for (std::size_t i = 0; i < 500; ++i) {
                auto pos = m::rnd::normal3(0., 1.);
                kernel->stateModel().addParticle({pos, idd});
            }

            readdy::scalar timeStep = 0.1;
            perform(kernel.get(), 1, timeStep);
        }
        SECTION("Should be rejected") {
            ctx.reactions().addEnzymatic("convAB", "C", "A", "B", 1e15, interactionDistance);
            ctx.reactions().addEnzymatic("convBA", "C", "B", "A", 1e-15, interactionDistance);
            ctx.particleTypes().add("D", 0.);
            ctx.potentials().addHarmonicRepulsion("B", "D", 10., 2.); // A and D interact, thus the initial state is unfavorable

            const auto idForward = ctx.reactions().idOf("convAB");
            const auto idBackward = ctx.reactions().idOf("convBA");

            auto countsObs = kernel->observe().reactionCounts(1,[&idForward, &idBackward](const readdy::model::observables::ReactionCounts::result_type &result) {
                const auto& counts = std::get<0>(result);
                if (counts.empty()) {
                    readdy::log::trace("reaction counts is empty, no reaction handler ran so far, skip test");
                    return;
                }
                REQUIRE(counts.at(idBackward) == 0);
                REQUIRE(counts.at(idForward) == 0);
            });
            auto countsConnection = kernel->connectObservable(countsObs.get());

            std::vector<std::string> typesToCount = {"A", "B", "C", "D"};
            auto numbersObs = kernel->observe().nParticles(1, typesToCount, [](const readdy::model::observables::NParticles::result_type &result) {
                REQUIRE(result[0] + result[1] == 1); // conservation of A + B
                REQUIRE(result[2] == 1); // conservation of C
                REQUIRE(result[3] == 500); // conservation of D
            });
            auto numbersConnection = kernel->connectObservable(numbersObs.get());

            const auto ida = ctx.particleTypes().idOf("A");
            const auto idc = ctx.particleTypes().idOf("C");
            const auto idd = ctx.particleTypes().idOf("D");
            kernel->stateModel().addParticle({{0., 0., 0.}, ida});
            kernel->stateModel().addParticle({{0., 0., 1.}, idc});

            // induce rejection by many repulsing particles in the vicinity
            for (std::size_t i = 0; i < 500; ++i) {
                auto pos = m::rnd::normal3(0., 1.);
                kernel->stateModel().addParticle({pos, idd});
            }

            readdy::scalar timeStep = 0.1;
            perform(kernel.get(), 1, timeStep);
        }
    }

    SECTION("Conservation of particles") {
        ctx.reactions().addFusion("fusion", "A", "B", "C", 10, reactionRadius);
        ctx.reactions().addFission("fission", "C", "A", "B", 2, reactionRadius);

        std::vector<std::string> typesToCount = {"A", "B", "C"};
        auto numbersObs = kernel->observe().nParticles(1, typesToCount, [](const readdy::model::observables::NParticles::result_type &result) {
            REQUIRE(result[0] + result[2] == 20); // conservation of A + C
            REQUIRE(result[1] + result[2] == 20); // conservation of B + C
            if (result[0] + result[2] != 1) {
                readdy::log::trace("A {} B {} C {}", result[0], result[1], result[2]);
            }
        });
        auto numbersConnection = kernel->connectObservable(numbersObs.get());

        const auto ida = ctx.particleTypes().idOf("A");
        const auto idb = ctx.particleTypes().idOf("B");
        const int n_particles = 20;
        std::vector<readdy::model::Particle> particlesA;
        std::vector<readdy::model::Particle> particlesB;
        for (auto i = 0; i < n_particles; ++i) {
            particlesA.emplace_back(readdy::model::rnd::normal3<readdy::scalar>(), ida);
            particlesB.emplace_back(readdy::model::rnd::normal3<readdy::scalar>(), idb);
        }
        kernel->stateModel().addParticles(particlesA);
        kernel->stateModel().addParticles(particlesB);

        readdy::scalar timeStep = 0.01;
        perform(kernel.get(), 100, timeStep, true);
    }

    SECTION("Mixed non-reversible and reversible") {
        ctx.reactions().addFusion("fusion", "A", "B", "C", 100, reactionRadius);
        ctx.reactions().addFission("fission", "C", "A", "B", 10, reactionRadius);
        ctx.reactions().add("nonRevConversion: A -> B", 10.);
        ctx.reactions().add("nonRevEnzymatic: B +(2) C -> A + C", 10.);
        // A + B + 2*C = const

        std::vector<std::string> typesToCount = {"A", "B", "C"};
        auto numbersObs = kernel->observe().nParticles(1, typesToCount, [](const readdy::model::observables::NParticles::result_type &result) {
            REQUIRE(result[0] + result[1] + 2*result[2] == 40); // conservation of A + B + 2C
        });
        auto numbersConnection = kernel->connectObservable(numbersObs.get());

        const auto ida = ctx.particleTypes().idOf("A");
        const auto idb = ctx.particleTypes().idOf("B");
        const int n_particles = 20;
        std::vector<readdy::model::Particle> particlesA;
        std::vector<readdy::model::Particle> particlesB;
        for (auto i = 0; i < n_particles; ++i) {
            particlesA.emplace_back(readdy::model::rnd::normal3<readdy::scalar>(), ida);
            particlesB.emplace_back(readdy::model::rnd::normal3<readdy::scalar>(), idb);
        }
        kernel->stateModel().addParticles(particlesA);
        kernel->stateModel().addParticles(particlesB);

        readdy::scalar timeStep = 0.1;
        perform(kernel.get(), 1000, timeStep, true);
    }
}
