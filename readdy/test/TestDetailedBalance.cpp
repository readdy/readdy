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
 * @copyright GPL-3
 */

#include <gtest/gtest.h>
#include <readdy/common/common.h>
#include <readdy/testing/KernelTest.h>
#include <readdy/testing/Utils.h>
#include <readdy/model/actions/Action.h>

namespace {

class TestDetailedBalanceWithKernels : public KernelTest {

};

namespace m = readdy::model;
namespace r = m::actions::reactions;

auto registerAbcSpecies = [](m::Context &ctx) {
    ctx.particleTypes().add("A", 1.);
    ctx.particleTypes().add("B", 1.);
    ctx.particleTypes().add("C", 1.);
};

TEST(TestDetailedBalance, ReversibleReactionConfigValuesFusionFission) {
    m::Context ctx;
    ctx.boxSize() = {21.544346900318832, 21.544346900318832, 21.544346900318832};
    ctx.kBT() = 1;
    registerAbcSpecies(ctx);
    ctx.potentials().addHarmonicRepulsion("A", "B", 2., 3.);
    ctx.reactions().add("fusion: A +(5) B -> C", 1.);
    ctx.reactions().add("fission: C -> A +(5) B", 2.);

    auto idFusion = ctx.reactions().idOf("fusion");
    auto idFission = ctx.reactions().idOf("fission");
    r::ReversibleReactionConfig conf(idFusion, idFission, ctx);
    EXPECT_NEAR(conf.totalVolume, 10000, 1e-4);
    EXPECT_FLOAT_EQ(conf.kbt, 1.);
    EXPECT_FLOAT_EQ(conf.lhsInteractionRadius, 3.);
    EXPECT_FLOAT_EQ(conf.reactionRadius, 5.);
    EXPECT_EQ(conf.numberLhsTypes, 2);
    EXPECT_EQ(conf.numberRhsTypes, 1);
    EXPECT_EQ(conf.rhsTypes[0], ctx.particleTypes().idOf("C"));
    EXPECT_EQ(conf.lhsTypes[0], ctx.particleTypes().idOf("A"));
    EXPECT_EQ(conf.lhsTypes[1], ctx.particleTypes().idOf("B"));
    EXPECT_NEAR(conf.lhsInteractionVolume, 113.09733, 1e-4);
    EXPECT_NEAR(conf.effectiveLhsInteractionVolume, 68.09910, 1e-4);
    EXPECT_NEAR(conf.effectiveLhsReactionVolume, 478.60057, 1e-4);
    EXPECT_NEAR(conf.equilibriumConstant, 41.60045, 1e-4);
    EXPECT_FLOAT_EQ(conf.macroBackwardRate, 2.);
    EXPECT_NEAR(conf.macroForwardRate, 480.76392, 1e-4);
    EXPECT_FLOAT_EQ(conf.acceptancePrefactor, 1.);
}

TEST(TestDetailedBalance, ReversibleReactionConfigValuesConversionConversion) {
    m::Context ctx;
    ctx.boxSize() = {10., 10., 10.};
    ctx.kBT() = 1;
    registerAbcSpecies(ctx);
    ctx.reactions().add("forward: A -> B", 1.);
    ctx.reactions().add("backward: B -> A", 2.);

    auto idForward = ctx.reactions().idOf("forward");
    auto idBackward = ctx.reactions().idOf("backward");
    r::ReversibleReactionConfig conf(idForward, idBackward, ctx);
    EXPECT_NEAR(conf.totalVolume, 1000, 1e-4);
    EXPECT_EQ(conf.numberLhsTypes, 1);
    EXPECT_EQ(conf.numberRhsTypes, 1);
    EXPECT_EQ(conf.lhsTypes[0], ctx.particleTypes().idOf("A"));
    EXPECT_EQ(conf.rhsTypes[0], ctx.particleTypes().idOf("B"));
    EXPECT_FLOAT_EQ(conf.lhsInteractionRadius, 0.);
    EXPECT_FLOAT_EQ(conf.lhsInteractionVolume, 0.);
    EXPECT_FLOAT_EQ(conf.rhsInteractionRadius, 0.);
    EXPECT_FLOAT_EQ(conf.rhsInteractionVolume, 0.);

    EXPECT_FLOAT_EQ(conf.effectiveLhsInteractionVolume, 0.);
    EXPECT_FLOAT_EQ(conf.effectiveLhsReactionVolume, 0.);
    EXPECT_FLOAT_EQ(conf.effectiveRhsInteractionVolume, 0.);
    EXPECT_FLOAT_EQ(conf.effectiveRhsReactionVolume, 0.);

    EXPECT_FLOAT_EQ(conf.acceptancePrefactor, 1.);

    EXPECT_FLOAT_EQ(conf.macroForwardRate, 1.);
    EXPECT_FLOAT_EQ(conf.macroBackwardRate, 2.);
    EXPECT_FLOAT_EQ(conf.equilibriumConstant, 2.);
}

TEST(TestDetailedBalance, ReversibleReactionConfigValuesEnzymaticEnzymatic) {
    m::Context ctx;
    ctx.boxSize() = {10., 10., 10.};
    ctx.kBT() = 1;
    registerAbcSpecies(ctx);
    ctx.potentials().addHarmonicRepulsion("A", "C", 2., 3.);
    ctx.potentials().addHarmonicRepulsion("B", "C", 2., 3.);
    ctx.reactions().add("forward: A +(5) C -> B + C", 1.);
    ctx.reactions().add("backward: B +(5) C -> A + C", 2.);

    auto idForward = ctx.reactions().idOf("forward");
    auto idBackward = ctx.reactions().idOf("backward");
    r::ReversibleReactionConfig conf(idForward, idBackward, ctx);
    EXPECT_NEAR(conf.totalVolume, 1000, 1e-4);
    EXPECT_EQ(conf.numberLhsTypes, 2);
    EXPECT_EQ(conf.numberRhsTypes, 2);
    EXPECT_EQ(conf.lhsTypes[0], ctx.particleTypes().idOf("A"));
    EXPECT_EQ(conf.rhsTypes[0], ctx.particleTypes().idOf("B"));
    EXPECT_EQ(conf.lhsTypes[1], ctx.particleTypes().idOf("C"));
    EXPECT_EQ(conf.rhsTypes[1], ctx.particleTypes().idOf("C"));
    EXPECT_FLOAT_EQ(conf.lhsInteractionRadius, 3.);
    EXPECT_FLOAT_EQ(conf.lhsInteractionVolume, 4. / 3. * readdy::util::numeric::pi<readdy::scalar>() * std::pow(3., 3));
    EXPECT_FLOAT_EQ(conf.rhsInteractionRadius, 3.);
    EXPECT_FLOAT_EQ(conf.rhsInteractionVolume, 4. / 3. * readdy::util::numeric::pi<readdy::scalar>() * std::pow(3., 3));

    EXPECT_NEAR(conf.effectiveLhsInteractionVolume, 68.09910, 1e-4);
    EXPECT_NEAR(conf.effectiveLhsReactionVolume, 478.60057, 1e-4);
    EXPECT_NEAR(conf.effectiveRhsInteractionVolume, 68.09910, 1e-4);
    EXPECT_NEAR(conf.effectiveRhsReactionVolume, 478.60057, 1e-4);

    EXPECT_NEAR(conf.acceptancePrefactor, 1., 1e-4) << "ratio of effective reaction volumina";

    EXPECT_NEAR(conf.macroForwardRate, 501.15150, 1e-4);
    EXPECT_NEAR(conf.macroBackwardRate, 1002.30300, 1e-4);
    EXPECT_FLOAT_EQ(conf.equilibriumConstant, 2.);

}

TEST(TestDetailedBalance, ReversibleReactionConfigFalseInput) {
    {
        // reaction radii do not match
        m::Context ctx;
        registerAbcSpecies(ctx);
        ctx.reactions().add("forward: A +(4) C -> B + C", 1.);
        ctx.reactions().add("backward: B +(5) C -> A + C", 2.);
        auto idForward = ctx.reactions().idOf("forward");
        auto idBackward = ctx.reactions().idOf("backward");
        EXPECT_THROW(r::ReversibleReactionConfig conf(idForward, idBackward, ctx), std::logic_error);
    }

    {
        // reactions are not reversible w.r.t. species
        m::Context ctx;
        registerAbcSpecies(ctx);
        ctx.reactions().add("forward: A -> B", 1.);
        ctx.reactions().add("backward: B -> C", 2.);
        auto idForward = ctx.reactions().idOf("forward");
        auto idBackward = ctx.reactions().idOf("backward");
        EXPECT_THROW(r::ReversibleReactionConfig conf(idForward, idBackward, ctx), std::logic_error);
    }

    {
        // reactions are not reversible w.r.t. species
        m::Context ctx;
        registerAbcSpecies(ctx);
        ctx.reactions().add("forward: A +(2) B -> C", 1.);
        ctx.reactions().add("backward: C -> A +(2) A", 2.);
        auto idForward = ctx.reactions().idOf("forward");
        auto idBackward = ctx.reactions().idOf("backward");
        EXPECT_THROW(r::ReversibleReactionConfig conf(idForward, idBackward, ctx), std::logic_error);
    }

    {
        // reactions are not reversible w.r.t. reaction radii
        m::Context ctx;
        registerAbcSpecies(ctx);
        ctx.reactions().add("forward: A +(2) B -> C", 1.);
        ctx.reactions().add("backward: C -> A +(3) B", 2.);
        auto idForward = ctx.reactions().idOf("forward");
        auto idBackward = ctx.reactions().idOf("backward");
        EXPECT_THROW(r::ReversibleReactionConfig conf(idForward, idBackward, ctx), std::logic_error);
    }
}

TEST(TestDetailedBalance, ReversibleReactionConfigValidCases) {
    {
        // FusionFission forward and backward might be swapped
        m::Context ctx;
        registerAbcSpecies(ctx);
        ctx.reactions().add("forward: C -> B +(2) A", 2.);
        ctx.reactions().add("backward: A +(2) B -> C", 1.);
        auto idForward = ctx.reactions().idOf("forward");
        auto idBackward = ctx.reactions().idOf("backward");
        EXPECT_NO_THROW(r::ReversibleReactionConfig conf(idForward, idBackward, ctx));
    }
    {
        // EnzymaticEnzymatic without potentials -> acceptancePrefactor = 1
        m::Context ctx;
        registerAbcSpecies(ctx);
        ctx.reactions().add("forward: A +(4) C -> B + C", 1.);
        ctx.reactions().add("backward: B +(4) C -> A + C", 2.);
        auto idForward = ctx.reactions().idOf("forward");
        auto idBackward = ctx.reactions().idOf("backward");
        r::ReversibleReactionConfig conf(idForward, idBackward, ctx);
        EXPECT_FLOAT_EQ(conf.acceptancePrefactor, 1.);
    }
}

TEST_P(TestDetailedBalanceWithKernels, SearchReversibleReactionsEnzymaticEnzymatic) {
    auto &ctx = kernel->context();
    registerAbcSpecies(ctx);
    ctx.reactions().add("forward: A +(4) C -> B + C", 1.);
    ctx.reactions().add("backward: B +(4) C -> A + C", 2.);
    auto &&reactions = kernel->actions().detailedBalance(1.);
    const auto &reversibleReactions = reactions->reversibleReactions();
    EXPECT_EQ(reversibleReactions.size(), 1);
    EXPECT_EQ(reversibleReactions.front()->reversibleType, r::ReversibleType::EnzymaticEnzymatic);
}

TEST_P(TestDetailedBalanceWithKernels, SearchReversibleReactionsFusionFission) {
    auto &ctx = kernel->context();
    registerAbcSpecies(ctx);
    ctx.reactions().add("fusion: A +(2) B -> C", 1.);
    ctx.reactions().add("fission: C -> A +(2) B", 2.);
    auto &&reactions = kernel->actions().detailedBalance(1.);
    const auto &reversibleReactions = reactions->reversibleReactions();
    EXPECT_EQ(reversibleReactions.size(), 1);
    EXPECT_EQ(reversibleReactions.front()->reversibleType, r::ReversibleType::FusionFission);
}

TEST_P(TestDetailedBalanceWithKernels, SearchReversibleReactionsFissionFusion) {
    auto &ctx = kernel->context();
    registerAbcSpecies(ctx);
    ctx.reactions().add("fission: C -> A +(2) B", 2.);
    ctx.reactions().add("fusion: B +(2) A -> C", 1.);
    auto &&reactions = kernel->actions().detailedBalance(1.);
    const auto &reversibleReactions = reactions->reversibleReactions();
    EXPECT_EQ(reversibleReactions.size(), 1);
    EXPECT_EQ(reversibleReactions.front()->reversibleType, r::ReversibleType::FusionFission);
}

TEST_P(TestDetailedBalanceWithKernels, SearchReversibleReactionsMichaelisMenten) {
    // one should be able to have a chain of reversible reactions like in Michaelis Menten type
    // reaction A + B <--> C <--> A + D
    // add a non-reversible, B --> D
    auto &ctx = kernel->context();
    registerAbcSpecies(ctx);
    ctx.particleTypes().add("D", 1.);
    ctx.reactions().add("fusion1: A +(2) B -> C", 2.);
    ctx.reactions().add("fission1: C -> A +(2) B", 1.);
    ctx.reactions().add("fission2: C -> A +(3) D", 3.);
    ctx.reactions().add("fusion2: A +(3) D -> C", 4.);
    ctx.reactions().add("shortcut: B -> D", 5.);
    
    auto &&reactions = kernel->actions().detailedBalance(1.);
    const auto &reversibleReactions = reactions->reversibleReactions();
    EXPECT_EQ(reversibleReactions.size(), 2);
    EXPECT_EQ(reversibleReactions[0]->reversibleType, r::ReversibleType::FusionFission);
    EXPECT_EQ(reversibleReactions[1]->reversibleType, r::ReversibleType::FusionFission);
    std::size_t index1, index2;
    if (reversibleReactions[0]->forwardName == "fusion1") {
        index1 = 0;
        index2 = 1;
    } else {
        index1 = 1;
        index2 = 0;
    }
    EXPECT_EQ(reversibleReactions[index1]->forwardName, "fusion1");
    EXPECT_EQ(reversibleReactions[index2]->forwardName, "fusion2");
    EXPECT_EQ(reversibleReactions[index1]->backwardName, "fission1");
    EXPECT_EQ(reversibleReactions[index2]->backwardName, "fission2");
}

TEST_P(TestDetailedBalanceWithKernels, SearchReversibleReactionsManyMixed) {
    // A + B <--> C, reversible
    // A + C --> B, non reversible
    // A + B <--> C + B, reversible, B is catalyst
    // B + A --> C + A, non-reversible, A is catalyst
    // A <--> B, reversible
    // C --> B, non-reversible
    auto &ctx = kernel->context();
    registerAbcSpecies(ctx);
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
    EXPECT_EQ(reversibleReactions.size(), 3);
    for (const auto &rev : reversibleReactions) {
        switch (rev->reversibleType) {
            case readdy::model::actions::reactions::FusionFission:
                EXPECT_EQ(rev->forwardName, "fusionRev");
                EXPECT_EQ(rev->backwardName, "fissionRev");
                break;
            case readdy::model::actions::reactions::ConversionConversion:
                EXPECT_TRUE((rev->forwardName == "convRev1" and rev->backwardName == "convRev2") or
                            (rev->forwardName == "convRev2" and rev->backwardName == "convRev1"));
                break;
            case readdy::model::actions::reactions::EnzymaticEnzymatic:
                EXPECT_TRUE((rev->forwardName == "enzRev1" and rev->backwardName == "enzRev2") or
                            (rev->forwardName == "enzRev2" and rev->backwardName == "enzRev1"));
                break;
            default:
                throw std::logic_error("Each type should be there");
        }
    }
}

TEST_P(TestDetailedBalanceWithKernels, SearchReversibleReactionsDuplicateReactions) {
    // A + B <--> C and another A + B --> C. It is ambiguous which reaction should be reversible --> throw
    auto &ctx = kernel->context();
    registerAbcSpecies(ctx);
    ctx.reactions().add("fusion: A +(2) B -> C", 2.);
    ctx.reactions().add("fission: C -> A +(2) B", 1.);
    ctx.reactions().add("anotherFusion: A +(2) B -> C", 2.);
    
    EXPECT_THROW({ auto &&reactions = kernel->actions().detailedBalance(1.); }, std::logic_error);
}

TEST_P(TestDetailedBalanceWithKernels, SearchReversibleReactionsWrongWeights) {
    auto &ctx = kernel->context();
    registerAbcSpecies(ctx);
    ctx.reactions().addFusion("fusion", "A", "B", "C", 1., 2., 0.5, 0.5);
    ctx.reactions().addFission("fission", "C", "A", "B", 1., 2., 0.3, 0.7);
    
    auto &&reactions = kernel->actions().detailedBalance(1.);
    const auto &reversibleReactions = reactions->reversibleReactions();
    EXPECT_EQ(reversibleReactions.size(), 0);
}


auto abcFusionFissionContext = [](readdy::model::Context &ctx, readdy::scalar rateOn, readdy::scalar rateOff) {
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
    ctx.reactions().addFission("fission", "C", "A", "B", rateOff, reactionRadius);
    ctx.reactions().addFusion("fusion", "A", "B", "C", rateOn, reactionRadius);
    
};



auto perform = [](readdy::model::Kernel *kernel, size_t nSteps, readdy::scalar timeStep, bool withIntegrator = false) {
    auto &&integrator = kernel->actions().eulerBDIntegrator(timeStep);
    auto &&forces = kernel->actions().calculateForces();
    using update_nl = readdy::model::actions::UpdateNeighborList;
    auto &&initNeighborList = kernel->actions().updateNeighborList(update_nl::Operation::init, 0);
    auto &&neighborList = kernel->actions().updateNeighborList(update_nl::Operation::update, 0);
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
};

TEST_P(TestDetailedBalanceWithKernels, FusionThatShouldBeRejected) {
    auto &ctx = kernel->context();
    abcFusionFissionContext(ctx, 1e15, 1e-15); // high on rate, small off rate
    ctx.particleTypes().add("D", 0.);
    ctx.potentials().addHarmonicRepulsion("C", "D", 10., 2.);

    const auto idfus = ctx.reactions().idOf("fusion");
    const auto idfis = ctx.reactions().idOf("fission");

    auto countsObs = kernel->observe().reactionCounts(1);
    countsObs->callback() = [&idfus, &idfis](const readdy::model::observables::ReactionCounts::result_type &result) {
        if (result.empty()) {
            readdy::log::trace("reaction counts is empty, no reaction handler ran so far, skip test");
            return;
        }
        EXPECT_EQ(result.at(idfus), 0) << "fusion shall not occur";
        EXPECT_EQ(result.at(idfis), 0) << "fission shall not occur";
    };
    auto countsConnection = kernel->connectObservable(countsObs.get());

    std::vector<std::string> typesToCount = {"A", "B", "C", "D"};
    auto numbersObs = kernel->observe().nParticles(1, typesToCount);
    numbersObs->callback() = [](const readdy::model::observables::NParticles::result_type &result) {
        EXPECT_EQ(result[0] + result[2], 1) << "conservation of A + C";
        EXPECT_EQ(result[1] + result[2], 1) << "conservation of B + C";
        EXPECT_EQ(result[3], 500) << "conservation of D";
        if (result[0] + result[2] != 1) {
            readdy::log::trace("A {} B {} C {}", result[0], result[1], result[2]);
        }
    };
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

TEST_P(TestDetailedBalanceWithKernels, FissionThatShouldBeRejected) {
    auto &ctx = kernel->context();
    abcFusionFissionContext(ctx, 1e-15, 1e15); // very high off rate
    ctx.particleTypes().add("D", 0.);
    ctx.potentials().addHarmonicRepulsion("A", "D", 10., 2.);
    ctx.potentials().addHarmonicRepulsion("B", "D", 10., 2.);

    const auto idfus = ctx.reactions().idOf("fusion");
    const auto idfis = ctx.reactions().idOf("fission");

    auto countsObs = kernel->observe().reactionCounts(1);
    countsObs->callback() = [&idfus, &idfis](const readdy::model::observables::ReactionCounts::result_type &result) {
        if (result.empty()) {
            readdy::log::trace("reaction counts is empty, no reaction handler ran so far, skip test");
            return;
        }
        EXPECT_EQ(result.at(idfis), 0) << "fission shall not occur";
        EXPECT_EQ(result.at(idfus), 0) << "fusion shall not occur";
    };
    auto countsConnection = kernel->connectObservable(countsObs.get());

    std::vector<std::string> typesToCount = {"A", "B", "C", "D"};
    auto numbersObs = kernel->observe().nParticles(1, typesToCount);
    numbersObs->callback() = [](const readdy::model::observables::NParticles::result_type &result) {
        EXPECT_EQ(result[0] + result[2], 1) << "conservation of A + C";
        EXPECT_EQ(result[1] + result[2], 1) << "conservation of B + C";
        EXPECT_EQ(result[3], 500) << "conservation of D";
        if (result[0] + result[2] != 1) {
            readdy::log::trace("A {} B {} C {}", result[0], result[1], result[2]);
        }
    };
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

TEST_P(TestDetailedBalanceWithKernels, ConservationOfParticles) {
    auto &ctx = kernel->context();
    abcFusionFissionContext(ctx, 10., 2.);

    std::vector<std::string> typesToCount = {"A", "B", "C"};
    auto numbersObs = kernel->observe().nParticles(1, typesToCount);
    numbersObs->callback() = [](const readdy::model::observables::NParticles::result_type &result) {
        EXPECT_EQ(result[0] + result[2], 20) << "conservation of A + C";
        EXPECT_EQ(result[1] + result[2], 20) << "conservation of B + C";
        if (result[0] + result[2] != 1) {
            readdy::log::trace("A {} B {} C {}", result[0], result[1], result[2]);
        }
    };
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

TEST_P(TestDetailedBalanceWithKernels, FusionThatShouldBeAccepted) {
    auto &ctx = kernel->context();
    abcFusionFissionContext(ctx, 1e15, 1e-15); // high on rate, small off rate
    ctx.particleTypes().add("D", 0.);
    ctx.potentials().addHarmonicRepulsion("A", "D", 10., 2.);
    ctx.potentials().addHarmonicRepulsion("B", "D", 10., 2.);

    const auto idfus = ctx.reactions().idOf("fusion");
    const auto idfis = ctx.reactions().idOf("fission");

    auto countsObs = kernel->observe().reactionCounts(1);
    countsObs->callback() = [&idfus, &idfis](const readdy::model::observables::ReactionCounts::result_type &result) {
        if (result.empty()) {
            readdy::log::trace("reaction counts is empty, no reaction handler ran so far, skip test");
            return;
        }
        EXPECT_EQ(result.at(idfus), 1) << "fusion must occur";
        EXPECT_EQ(result.at(idfis), 0) << "fission shall not occur";
    };
    auto countsConnection = kernel->connectObservable(countsObs.get());

    std::vector<std::string> typesToCount = {"A", "B", "C", "D"};
    auto numbersObs = kernel->observe().nParticles(1, typesToCount);
    numbersObs->callback() = [](const readdy::model::observables::NParticles::result_type &result) {
        EXPECT_EQ(result[0] + result[2], 1) << "conservation of A + C";
        EXPECT_EQ(result[1] + result[2], 1) << "conservation of B + C";
        EXPECT_EQ(result[3], 500) << "conservation of D";
        if (result[0] + result[2] != 1) {
            readdy::log::trace("A {} B {} C {}", result[0], result[1], result[2]);
        }
    };
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

TEST_P(TestDetailedBalanceWithKernels, FissionThatShouldBeAccepted) {
    auto &ctx = kernel->context();
    abcFusionFissionContext(ctx, 1e-15, 1e15); // very high off rate
    ctx.particleTypes().add("D", 0.);
    ctx.potentials().addHarmonicRepulsion("C", "D", 10., 2.);

    const auto idfus = ctx.reactions().idOf("fusion");
    const auto idfis = ctx.reactions().idOf("fission");

    auto countsObs = kernel->observe().reactionCounts(1);
    countsObs->callback() = [&idfus, &idfis](const readdy::model::observables::ReactionCounts::result_type &result) {
        if (result.empty()) {
            readdy::log::trace("reaction counts is empty, no reaction handler ran so far, skip test");
            return;
        }
        EXPECT_EQ(result.at(idfis), 1) << "fission must occur";
        EXPECT_EQ(result.at(idfus), 0) << "fusion shall not occur";
    };
    auto countsConnection = kernel->connectObservable(countsObs.get());

    std::vector<std::string> typesToCount = {"A", "B", "C", "D"};
    auto numbersObs = kernel->observe().nParticles(1, typesToCount);
    numbersObs->callback() = [](const readdy::model::observables::NParticles::result_type &result) {
        EXPECT_EQ(result[0] + result[2], 1) << "conservation of A + C";
        EXPECT_EQ(result[1] + result[2], 1) << "conservation of B + C";
        EXPECT_EQ(result[3], 500) << "conservation of D";
        if (result[0] + result[2] != 1) {
            readdy::log::trace("A {} B {} C {}", result[0], result[1], result[2]);
        }
    };
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

TEST_P(TestDetailedBalanceWithKernels, MixedNonReversibleAndReversible) {
    auto &ctx = kernel->context();
    abcFusionFissionContext(ctx, 100., 10.);
    ctx.reactions().add("nonRevConversion: A -> B", 10.);
    ctx.reactions().add("nonRevEnzymatic: B +(2) C -> A + C", 10.);
    // A + B + 2*C = const

    std::vector<std::string> typesToCount = {"A", "B", "C"};
    auto numbersObs = kernel->observe().nParticles(1, typesToCount);
    numbersObs->callback() = [](const readdy::model::observables::NParticles::result_type &result) {
        EXPECT_EQ(result[0] + result[1] + 2*result[2], 40) << "conservation of A + B + 2C";
    };
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

auto abConversionContext = [](readdy::model::Context &ctx, readdy::scalar rateOn, readdy::scalar rateOff) {
    ctx.boxSize() = {{12, 12, 12}};
    ctx.particleTypes().add("A", 1.);
    ctx.particleTypes().add("B", 1.);
    readdy::scalar interactionDistance = 2.;
    ctx.potentials().addHarmonicRepulsion("A", "B", 10., interactionDistance);
    ctx.potentials().addHarmonicRepulsion("B", "B", 10., interactionDistance);
    ctx.potentials().addHarmonicRepulsion("A", "A", 10., interactionDistance);
    ctx.reactions().addConversion("convAB", "A", "B", rateOn);
    ctx.reactions().addConversion("convBA", "B", "A", rateOff);
    
};

TEST_P(TestDetailedBalanceWithKernels, ConversionThatShouldBeAccepted) {
    auto &ctx = kernel->context();
    abConversionContext(ctx, 1e15, 1e-15);
    ctx.particleTypes().add("D", 0.);
    ctx.potentials().addHarmonicRepulsion("A", "D", 10., 2.); // A and D interact, thus the initial state is unfavorable

    const auto idForward = ctx.reactions().idOf("convAB");
    const auto idBackward = ctx.reactions().idOf("convBA");

    auto countsObs = kernel->observe().reactionCounts(1);
    countsObs->callback() = [&idForward, &idBackward](const readdy::model::observables::ReactionCounts::result_type &result) {
        if (result.empty()) {
            readdy::log::trace("reaction counts is empty, no reaction handler ran so far, skip test");
            return;
        }
        EXPECT_EQ(result.at(idBackward), 0);
        EXPECT_EQ(result.at(idForward), 1);
    };
    auto countsConnection = kernel->connectObservable(countsObs.get());

    std::vector<std::string> typesToCount = {"A", "B", "D"};
    auto numbersObs = kernel->observe().nParticles(1, typesToCount);
    numbersObs->callback() = [](const readdy::model::observables::NParticles::result_type &result) {
        EXPECT_EQ(result[0] + result[1], 1) << "conservation of A + B";
        EXPECT_EQ(result[2], 500) << "conservation of D";
    };
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

TEST_P(TestDetailedBalanceWithKernels, ConversionThatShouldBeRejected) {
    auto &ctx = kernel->context();
    abConversionContext(ctx, 1e15, 1e-15);
    ctx.particleTypes().add("D", 0.);
    ctx.potentials().addHarmonicRepulsion("B", "D", 10., 2.);

    const auto idForward = ctx.reactions().idOf("convAB");
    const auto idBackward = ctx.reactions().idOf("convBA");

    auto countsObs = kernel->observe().reactionCounts(1);
    countsObs->callback() = [&idForward, &idBackward](const readdy::model::observables::ReactionCounts::result_type &result) {
        if (result.empty()) {
            readdy::log::trace("reaction counts is empty, no reaction handler ran so far, skip test");
            return;
        }
        EXPECT_EQ(result.at(idBackward), 0);
        EXPECT_EQ(result.at(idForward), 0);
    };
    auto countsConnection = kernel->connectObservable(countsObs.get());

    std::vector<std::string> typesToCount = {"A", "B", "D"};
    auto numbersObs = kernel->observe().nParticles(1, typesToCount);
    numbersObs->callback() = [](const readdy::model::observables::NParticles::result_type &result) {
        EXPECT_EQ(result[0] + result[1], 1) << "conservation of A + B";
        EXPECT_EQ(result[2], 500) << "conservation of D";
    };
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

auto abEnzymaticContext = [](readdy::model::Context &ctx, readdy::scalar rateOn, readdy::scalar rateOff) {
    ctx.boxSize() = {{12, 12, 12}};
    ctx.particleTypes().add("A", 1.);
    ctx.particleTypes().add("B", 1.);
    ctx.particleTypes().add("C", 1.);
    readdy::scalar interactionDistance = 2.;
    ctx.potentials().addHarmonicRepulsion("A", "B", 10., interactionDistance);
    ctx.potentials().addHarmonicRepulsion("B", "B", 10., interactionDistance);
    ctx.potentials().addHarmonicRepulsion("A", "A", 10., interactionDistance);
    ctx.reactions().addEnzymatic("convAB", "C", "A", "B", rateOn, interactionDistance);
    ctx.reactions().addEnzymatic("convBA", "C", "B", "A", rateOff, interactionDistance);
    
};

TEST_P(TestDetailedBalanceWithKernels, EnzymaticThatShouldBeAccepted) {
    auto &ctx = kernel->context();
    abEnzymaticContext(ctx, 1e15, 1e-15);
    ctx.particleTypes().add("D", 0.);
    ctx.potentials().addHarmonicRepulsion("A", "D", 10., 2.); // A and D interact, thus the initial state is unfavorable

    const auto idForward = ctx.reactions().idOf("convAB");
    const auto idBackward = ctx.reactions().idOf("convBA");

    auto countsObs = kernel->observe().reactionCounts(1);
    countsObs->callback() = [&idForward, &idBackward](const readdy::model::observables::ReactionCounts::result_type &result) {
        if (result.empty()) {
            readdy::log::trace("reaction counts is empty, no reaction handler ran so far, skip test");
            return;
        }
        EXPECT_EQ(result.at(idBackward), 0);
        EXPECT_EQ(result.at(idForward), 1);
    };
    auto countsConnection = kernel->connectObservable(countsObs.get());

    std::vector<std::string> typesToCount = {"A", "B", "C", "D"};
    auto numbersObs = kernel->observe().nParticles(1, typesToCount);
    numbersObs->callback() = [](const readdy::model::observables::NParticles::result_type &result) {
        EXPECT_EQ(result[0] + result[1], 1) << "conservation of A + B";
        EXPECT_EQ(result[2], 1) << "conservation of C";
        EXPECT_EQ(result[3], 500) << "conservation of D";
    };
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

TEST_P(TestDetailedBalanceWithKernels, EnzymaticThatShouldBeRejected) {
    auto &ctx = kernel->context();
    abEnzymaticContext(ctx, 1e15, 1e-15);
    ctx.particleTypes().add("D", 0.);
    ctx.potentials().addHarmonicRepulsion("B", "D", 10., 2.); // A and D interact, thus the initial state is unfavorable

    const auto idForward = ctx.reactions().idOf("convAB");
    const auto idBackward = ctx.reactions().idOf("convBA");

    auto countsObs = kernel->observe().reactionCounts(1);
    countsObs->callback() = [&idForward, &idBackward](const readdy::model::observables::ReactionCounts::result_type &result) {
        if (result.empty()) {
            readdy::log::trace("reaction counts is empty, no reaction handler ran so far, skip test");
            return;
        }
        EXPECT_EQ(result.at(idBackward), 0);
        EXPECT_EQ(result.at(idForward), 0);
    };
    auto countsConnection = kernel->connectObservable(countsObs.get());

    std::vector<std::string> typesToCount = {"A", "B", "C", "D"};
    auto numbersObs = kernel->observe().nParticles(1, typesToCount);
    numbersObs->callback() = [](const readdy::model::observables::NParticles::result_type &result) {
        EXPECT_EQ(result[0] + result[1], 1) << "conservation of A + B";
        EXPECT_EQ(result[2], 1) << "conservation of C";
        EXPECT_EQ(result[3], 500) << "conservation of D";
    };
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

std::vector<std::string> singlecpu = {{"SingleCPU"}};

INSTANTIATE_TEST_CASE_P(TestDetailedBalanceWithKernels, TestDetailedBalanceWithKernels,
//::testing::ValuesIn(readdy::testing::getKernelsToTest()));
                        ::testing::ValuesIn(singlecpu));

}
