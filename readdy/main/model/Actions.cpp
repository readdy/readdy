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
 * @file Actions.cpp
 * @brief Core library implementation of Actions
 * @author clonker
 * @author chrisfroe
 * @date 17.11.16
 */

#include <readdy/model/actions/Actions.h>
#include <readdy/model/Kernel.h>
#include <readdy/common/integration.h>
#include <readdy/common/numeric.h>
#include <readdy/model/actions/DetailedBalance.h>

namespace readdy {
namespace model {
namespace actions {


UpdateNeighborList::UpdateNeighborList(UpdateNeighborList::Operation operation, scalar skinSize)
        : operation(operation), skinSize(skinSize) {
}

EulerBDIntegrator::EulerBDIntegrator(scalar timeStep) : TimeStepDependentAction(timeStep) {}

reactions::UncontrolledApproximation::UncontrolledApproximation(scalar timeStep) : TimeStepDependentAction(timeStep) {}

reactions::Gillespie::Gillespie(scalar timeStep) : TimeStepDependentAction(timeStep) {}

reactions::ReversibleReactionConfig::ReversibleReactionConfig(ReactionId forwardId, ReactionId backwardId,
                                                              const Context &ctx)
        : forwardId(forwardId), backwardId(backwardId) {

    forwardReaction = ctx.reactions().byId(forwardId);
    backwardReaction = ctx.reactions().byId(backwardId);

    // species
    numberLhsTypes = forwardReaction->nEducts();
    lhsTypes = forwardReaction->educts();
    numberRhsTypes = backwardReaction->nEducts();
    rhsTypes = backwardReaction->educts();
    if (numberLhsTypes == 2 and numberRhsTypes == 1) {
        reversibleType = ReversibleType::FusionFission;
        lhsNames = {ctx.particleTypes().nameOf(lhsTypes[0]), ctx.particleTypes().nameOf(lhsTypes[1])};
        rhsNames = {ctx.particleTypes().nameOf(rhsTypes[0]), ""};
    } else if (numberLhsTypes == 1 and numberRhsTypes == 2) {
        reversibleType = ReversibleType::FusionFission;
        // in FusionFission the lhs is assumed to have two educts -> swap forward and backward and redetermine types
        std::swap(forwardReaction, backwardReaction);
        numberLhsTypes = forwardReaction->nEducts();
        lhsTypes = forwardReaction->educts();
        numberRhsTypes = backwardReaction->nEducts();
        rhsTypes = backwardReaction->educts();
        lhsNames = {ctx.particleTypes().nameOf(lhsTypes[0]), ctx.particleTypes().nameOf(lhsTypes[1])};
        rhsNames = {ctx.particleTypes().nameOf(rhsTypes[0]), ""};
    } else if (numberLhsTypes == 1 and numberRhsTypes == 1) {
        reversibleType = ReversibleType::ConversionConversion;
        lhsNames = {ctx.particleTypes().nameOf(lhsTypes[0]), ""};
        rhsNames = {ctx.particleTypes().nameOf(rhsTypes[0]), ""};
    } else if (numberLhsTypes == 2 and numberRhsTypes == 2) {
        reversibleType = ReversibleType::EnzymaticEnzymatic;
        lhsNames = {ctx.particleTypes().nameOf(lhsTypes[0]), ctx.particleTypes().nameOf(lhsTypes[1])};
        rhsNames = {ctx.particleTypes().nameOf(rhsTypes[0]), ctx.particleTypes().nameOf(rhsTypes[1])};
    } else {
        throw std::logic_error(fmt::format("Type of reversible reaction cannot be determined, method: {} file: {}",
                                           "ReversibleReactionConfig::ReversibleReactionConfig", "Actions.cpp"));
    }

    // check for reversibility in particle species
    if (not(lhsTypes == backwardReaction->products()) or not(rhsTypes == forwardReaction->products())) {
        if (reversibleType == ReversibleType::FusionFission and
            (lhsTypes[0] == backwardReaction->products()[1] and lhsTypes[1] == backwardReaction->products()[0]) and
                (forwardReaction->products()[0] == backwardReaction->educts()[0])) {
            // FusionFission might have swapped order of lhs types forward and backward, which is fine
        } else {
            throw std::logic_error(fmt::format(
                    "The given reactions are not reversible with respect to particle species, method: {} file: {}",
                    "ReversibleReactionConfig::ReversibleReactionConfig", "Actions.cpp"));
        }
    }

    // check for reversibility in weights for FusionFission
    if (reversibleType == ReversibleType::FusionFission) {
        if (forwardReaction->educts()[0] == backwardReaction->products()[0] and
            forwardReaction->educts()[1] == backwardReaction->products()[1]) {
            // conventional order, forward weights must be backward weights
            if (not (forwardReaction->weight1()==backwardReaction->weight1() and
                     forwardReaction->weight2()==backwardReaction->weight2())) {
                throw std::logic_error(fmt::format(
                        "The given reactions are not reversible with respect to weights, method: {} file: {}",
                        "ReversibleReactionConfig::ReversibleReactionConfig", "Actions.cpp"));
            }
        } else if (forwardReaction->educts()[0] == backwardReaction->products()[1] and
                   forwardReaction->educts()[1] == backwardReaction->products()[0]) {
            // flipped order, forward weights must be flipped backward weights
            if (not (forwardReaction->weight1()==backwardReaction->weight2() and
                     forwardReaction->weight2()==backwardReaction->weight1())) {
                throw std::logic_error(fmt::format(
                        "The given reactions are not reversible with respect to weights, method: {} file: {}",
                        "ReversibleReactionConfig::ReversibleReactionConfig", "Actions.cpp"));
            }
        }
    }

    forwardName = forwardReaction->name();
    backwardName = backwardReaction->name();

    microForwardRate = forwardReaction->rate();
    microBackwardRate = backwardReaction->rate();

    // reaction radius
    switch (reversibleType) {
        case FusionFission:
            if (forwardReaction->eductDistance() != backwardReaction->productDistance()) {
                throw std::logic_error(fmt::format("Reaction radii forward and backward must be equal,"
                                                   "were {} and {}, method: {} file: {}",
                                                   forwardReaction->eductDistance(),
                                                   backwardReaction->productDistance(),
                                                   "ReversibleReactionConfig::ReversibleReactionConfig",
                                                   "Actions.cpp"));
            }
            reactionRadius = forwardReaction->eductDistance();
            break;
        case ConversionConversion:
            // no action required
            break;
        case EnzymaticEnzymatic:
            // Note that Enzymatics only define the educt distance
            if (forwardReaction->eductDistance() != backwardReaction->eductDistance()) {
                throw std::logic_error(fmt::format("Reaction radii forward and backward must be equal,"
                                                   "were {} and {}, method: {} file: {}",
                                                   forwardReaction->eductDistance(),
                                                   backwardReaction->eductDistance(),
                                                   "ReversibleReactionConfig::ReversibleReactionConfig",
                                                   "Actions.cpp"));
            }
            reactionRadius = forwardReaction->eductDistance();
            break;
        default:
            throw std::logic_error(fmt::format("Unknown type of reversible reaction, method: {} file: {}",
                                               "ReversibleReactionConfig::ReversibleReactionConfig", "Actions.cpp"));
    }

    // potentials
    switch(reversibleType) {
        case FusionFission:
            lhsPotentials = ctx.potentials().potentialsOf(lhsTypes[0], lhsTypes[1]);
            break;
        case ConversionConversion:
            // no action required
            break;
        case EnzymaticEnzymatic:
            lhsPotentials = ctx.potentials().potentialsOf(lhsTypes[0], lhsTypes[1]);
            rhsPotentials = ctx.potentials().potentialsOf(rhsTypes[0], rhsTypes[1]);
            break;
        default:
            throw std::logic_error(fmt::format("Unknown type of reversible reaction, method: {} file: {}",
                                               "ReversibleReactionConfig::ReversibleReactionConfig", "Actions.cpp"));
    }

    // interaction radius and volume
    if (lhsPotentials.empty()) {
        lhsInteractionRadius = 0.;
        lhsInteractionVolume = 0.;
    } else {
        auto maxLhsCutoffPotential = std::max_element(lhsPotentials.begin(), lhsPotentials.end(),
                                                      [](const pot p1, const pot p2) {
                                                          return p1->getCutoffRadius() < p2->getCutoffRadius();
                                                      });
        lhsInteractionRadius = (*maxLhsCutoffPotential)->getCutoffRadius();
        lhsInteractionVolume = 4./3. * util::numeric::pi<scalar>() * std::pow(lhsInteractionRadius, 3);
    }

    if (rhsPotentials.empty()) {
        rhsInteractionRadius = 0.;
        rhsInteractionVolume = 0.;
    } else {
        auto maxRhsCutoffPotential = std::max_element(rhsPotentials.begin(), rhsPotentials.end(),
                                                      [](const pot p1, const pot p2) {
                                                          return p1->getCutoffRadius() < p2->getCutoffRadius();
                                                      });
        rhsInteractionRadius = (*maxRhsCutoffPotential)->getCutoffRadius();
        rhsInteractionVolume = 4./3. * util::numeric::pi<scalar>() * std::pow(rhsInteractionRadius, 3);
    }

    // effective interaction and reaction volumina
    const auto integrand = [](const scalar radius, std::vector<pot> potentials, const scalar kbt_) -> scalar {
        scalar energy = 0.;
        for (const auto &pot : potentials) {
            const Vec3 difference(radius, 0., 0.);
            energy += pot->calculateEnergy(difference);
        }
        const scalar result =
                4. * readdy::util::numeric::pi<scalar>() * radius * radius * std::exp(-1. * energy / kbt_);
        return result;
    };
    kbt = ctx.kBT();
    const auto lhsIntegrand = std::bind(integrand, std::placeholders::_1, lhsPotentials, kbt);
    const auto rhsIntegrand = std::bind(integrand, std::placeholders::_1, rhsPotentials, kbt);

    const scalar desiredRelativeError = 1e-12;
    const std::size_t maxiter = 1000;

    // effective lhs interaction volume
    {
        const auto result = util::integration::integrateAdaptive(lhsIntegrand, 0., lhsInteractionRadius,
                                                                 desiredRelativeError, maxiter);
        const auto achievedRelativeError = result.second / result.first;
        if (achievedRelativeError > desiredRelativeError) {
            readdy::log::warn("Integration of potentials yielded larger error ({}) than desired ({}).",
                              achievedRelativeError, desiredRelativeError);
        }
        effectiveLhsInteractionVolume = result.first;
    }

    // effective rhs interaction volume
    {
        const auto result = util::integration::integrateAdaptive(rhsIntegrand, 0., rhsInteractionRadius,
                                                                 desiredRelativeError, maxiter);
        const auto achievedRelativeError = result.second / result.first;
        if (achievedRelativeError> desiredRelativeError) {
            readdy::log::warn("Integration of potentials yielded larger error ({}) than desired ({}).",
                              achievedRelativeError, desiredRelativeError);
        }
        effectiveRhsInteractionVolume = result.first;
    }

    // effective lhs reaction volume
    {
        const auto result = util::integration::integrateAdaptive(lhsIntegrand, 0., reactionRadius, desiredRelativeError,
                                                                 maxiter);
        const auto achievedRelativeError = result.second / result.first;
        if (achievedRelativeError> desiredRelativeError) {
            readdy::log::warn("Integration of potentials yielded larger error ({}) than desired ({}).",
                              achievedRelativeError, desiredRelativeError);
        }
        effectiveLhsReactionVolume = result.first;
    }

    // effective rhs reaction volume
    {
        const auto result = util::integration::integrateAdaptive(rhsIntegrand, 0., reactionRadius, desiredRelativeError,
                                                                 maxiter);
        const auto achievedRelativeError = result.second / result.first;
        if (achievedRelativeError> desiredRelativeError) {
            readdy::log::warn("Integration of potentials yielded larger error ({}) than desired ({}).",
                              achievedRelativeError, desiredRelativeError);
        }
        effectiveRhsReactionVolume = result.first;
    }

    // FusionFission needs a cumulative distribution
    if (reversibleType == ReversibleType::FusionFission) {
        const size_t nPoints = 1000;
        fissionRadii.reserve(nPoints);
        const scalar dr = reactionRadius / static_cast<scalar>(nPoints);
        scalar r = 0.;
        for (size_t i = 0; i < nPoints + 1; ++i) {
            fissionRadii.push_back(r);
            r += dr;
        }
        for (size_t i = 0; i < fissionRadii.size(); ++i) {
            const auto result = util::integration::integrateAdaptive(lhsIntegrand, fissionRadii.front(),
                                                                     fissionRadii[i], desiredRelativeError, maxiter);
            const auto cumvalue = result.first;
            const auto achievedRelativeError = result.second /  result.first;
            if (achievedRelativeError > desiredRelativeError) {
                readdy::log::warn("Integration of lhs potentials yielded larger error ({}) than desired ({}).",
                                  achievedRelativeError, desiredRelativeError);
            }
            cumulativeFissionProb.push_back(cumvalue);
        }
        if (cumulativeFissionProb.back() >= 1e-16) {
            // normalize fission probability to 1
            for (auto &p : cumulativeFissionProb) {
                p = p / cumulativeFissionProb.back();
            }
        } else {
            // in case of small values, construct a cumulative that is zero everywhere and 1 in the last entry 
            std::for_each(cumulativeFissionProb.begin(), cumulativeFissionProb.end() - 1,
                          [](auto &entry) { entry = 0.; });
            cumulativeFissionProb.back() = 1.;
        }
    }

    // acceptance prefactor
    switch (reversibleType) {
        case FusionFission:
            acceptancePrefactor = 1.;
            break;
        case ConversionConversion:
            acceptancePrefactor = 1.;
            break;
        case EnzymaticEnzymatic:
            acceptancePrefactor = effectiveLhsReactionVolume / effectiveRhsReactionVolume;
            break;
        default:
            throw std::logic_error(fmt::format("Unknown type of reversible reaction, method: {} file: {}",
                                               "ReversibleReactionConfig::ReversibleReactionConfig", "Actions.cpp"));
    }

    totalVolume = ctx.boxVolume();

    // inferred macro rates
    switch (reversibleType) {
        case FusionFission:
            // In this case the equ. constant is K_d V, where K_d = k_off / k_on is the dissociation constant.
            // We multiply it by the volume to have a unitless expression
            equilibriumConstant = (microBackwardRate / microForwardRate) *
                                  (totalVolume - lhsInteractionVolume + effectiveLhsInteractionVolume) /
                                  effectiveLhsReactionVolume;
            macroBackwardRate = microBackwardRate;
            macroForwardRate = macroBackwardRate * totalVolume / equilibriumConstant;
            break;
        case ConversionConversion:
            macroForwardRate = microForwardRate;
            macroBackwardRate = microBackwardRate;
            equilibriumConstant = macroBackwardRate / macroForwardRate;
            break;
        case EnzymaticEnzymatic:
            macroForwardRate = microForwardRate * totalVolume * effectiveLhsReactionVolume /
                               (totalVolume - lhsInteractionVolume + effectiveLhsInteractionVolume);
            macroBackwardRate = microBackwardRate * totalVolume * effectiveRhsReactionVolume /
                                (totalVolume - rhsInteractionVolume + effectiveRhsInteractionVolume);
            equilibriumConstant = macroBackwardRate / macroForwardRate;
            break;
        default:
            throw std::logic_error(fmt::format("Unknown type of reversible reaction, method: {} file: {}",
                                               "ReversibleReactionConfig::ReversibleReactionConfig", "Actions.cpp"));
    }
}

std::string reactions::ReversibleReactionConfig::describe() const {
    std::stringstream description;

    description << fmt::format("ReversibleReactionConfig({}", readdy::util::str::newline);
    description << fmt::format("reversibleReactionType {} {}", reversibleType, readdy::util::str::newline);

    if (numberLhsTypes == 1) {
        description
                << fmt::format("lhsTypes {} {}", lhsNames[0], readdy::util::str::newline);
    } else if (numberLhsTypes == 2) {
        description << fmt::format("lhsTypes {} and {} {}", lhsNames[0], lhsNames[1], readdy::util::str::newline);
    }
    if (numberRhsTypes == 1) {
        description
                << fmt::format("rhsTypes {} {}", rhsNames[0], readdy::util::str::newline);
    } else if (numberRhsTypes == 2) {
        description << fmt::format("rhsTypes {} and {} {}", rhsNames[0], rhsNames[1], readdy::util::str::newline);
    }

    description << fmt::format("forwardReaction \"{}\"{}{}{}", forwardName,
                               readdy::util::str::newline, *forwardReaction, readdy::util::str::newline);
    description << fmt::format("backwardReaction \"{}\"{}{}{}", backwardName,
                               readdy::util::str::newline, *backwardReaction, readdy::util::str::newline);
    description << fmt::format("forwardId {} {}", forwardId, readdy::util::str::newline);
    description << fmt::format("backwardId {} {}", backwardId, readdy::util::str::newline);

    description << fmt::format("microForwardRate {} {}", microForwardRate, readdy::util::str::newline);
    description << fmt::format("microBackwardRate {} {}", microBackwardRate, readdy::util::str::newline);
    description << fmt::format("reactionRadius {} {}", reactionRadius, readdy::util::str::newline);

    description << fmt::format("lhsPotentials{}", readdy::util::str::newline);
    for (const auto &p : lhsPotentials) {
        description << fmt::format("    {}{}", p->describe(), readdy::util::str::newline);
    }
    description << fmt::format("rhsPotentials{}", readdy::util::str::newline);
    for (const auto &p : rhsPotentials) {
        description << fmt::format("    {}{}", p->describe(), readdy::util::str::newline);
    }
    
    description << fmt::format("lhsInteractionRadius {} {}", lhsInteractionRadius, readdy::util::str::newline);
    description << fmt::format("lhsInteractionVolume {} {}", lhsInteractionVolume, readdy::util::str::newline);
    description << fmt::format("effectiveLhsInteractionVolume {} {}", effectiveLhsInteractionVolume, readdy::util::str::newline);
    description << fmt::format("effectiveLhsReactionVolume {} {}", effectiveLhsReactionVolume, readdy::util::str::newline);

    description << fmt::format("rhsInteractionRadius {} {}", rhsInteractionRadius, readdy::util::str::newline);
    description << fmt::format("rhsInteractionVolume {} {}", rhsInteractionVolume, readdy::util::str::newline);
    description << fmt::format("effectiveRhsInteractionVolume {} {}", effectiveRhsInteractionVolume, readdy::util::str::newline);
    description << fmt::format("effectiveRhsReactionVolume {} {}", effectiveRhsReactionVolume, readdy::util::str::newline);
    
    description << fmt::format("totalVolume {} {}", totalVolume, readdy::util::str::newline);
    description << fmt::format("kbt {} {}", kbt, readdy::util::str::newline);
    description << fmt::format("acceptancePrefactor {} {}", acceptancePrefactor, readdy::util::str::newline);

    description << fmt::format("inferred macro rates {}", readdy::util::str::newline);
    description << fmt::format("    equilibriumConstant {} {}", equilibriumConstant, readdy::util::str::newline);
    description << fmt::format("    macroForwardRate {} {}", macroForwardRate, readdy::util::str::newline);
    description << fmt::format("    macroBackwardRate {} {}", macroBackwardRate, readdy::util::str::newline);
    description << fmt::format("){}",  readdy::util::str::newline);
    return description.str();
}

reactions::DetailedBalance::DetailedBalance(scalar timeStep) : TimeStepDependentAction(timeStep) {}

void reactions::DetailedBalance::searchReversibleReactions(const Context &ctx) {
    _reversibleReactionsMap.clear();
    _reversibleReactionsContainer.clear();

    const auto &order1Flat = ctx.reactions().order1Flat();
    const auto &order2Flat = ctx.reactions().order2Flat();

    // look for FusionFission
    for (const auto &reactionO2 : order2Flat) {
        for (const auto &reactionO1 : order1Flat) {
            if (reactionO2->nProducts() == 1 and reactionO1->nProducts() == 2 and
                reactionO1->productDistance() == reactionO2->eductDistance()) {
                // also consider as reversible if - for one reaction - the educts are swapped
                // for FusionFission also the weights must fit
                if ((reactionO1->educts()[0] == reactionO2->products()[0]) and (
                        (
                            (reactionO2->educts()[0] == reactionO1->products()[0]) and
                            (reactionO2->educts()[1] == reactionO1->products()[1]) and
                            (reactionO1->weight1() == reactionO2->weight1()) and
                            (reactionO1->weight2() == reactionO2->weight2())
                        ) or (
                            (reactionO2->educts()[0] == reactionO1->products()[1]) and
                            (reactionO2->educts()[0] == reactionO1->products()[1]) and
                                (reactionO1->weight1() == reactionO2->weight2()) and
                                (reactionO1->weight2() == reactionO2->weight1())
                        )
                )) {
                    _reversibleReactionsContainer.emplace_back(
                            std::make_shared<ReversibleReactionConfig>(reactionO2->id(), reactionO1->id(), ctx)
                    );
                }
            }
        }
    }

    // look for ConversionConversion
    for (std::size_t i = 0; i < order1Flat.size(); ++i) {
        for (std::size_t j = i + 1; j < order1Flat.size(); ++j) {
            const auto &ri = order1Flat[i];
            const auto &rj = order1Flat[j];
            if (ri->nProducts() == 1 and rj->nProducts() == 1) {
                if ((ri->educts()[0] == rj->products()[0]) and (ri->products()[0] == rj->educts()[0])) {
                    _reversibleReactionsContainer.emplace_back(
                            std::make_shared<ReversibleReactionConfig>(ri->id(), rj->id(), ctx)
                    );
                }
            }
        }
    }

    // look for EnzymaticEnzymatic
    for (std::size_t i = 0; i < order2Flat.size(); ++i) {
        for (std::size_t j = i + 1; j < order2Flat.size(); ++j) {
            const auto &ri = order2Flat[i];
            const auto &rj = order2Flat[j];
            if (ri->nProducts() == 2 and rj->nProducts() == 2 and ri->eductDistance()==rj->eductDistance()) {
                if ((ri->educts() == rj->products()) and (ri->products() == rj->educts())) {
                    _reversibleReactionsContainer.emplace_back(
                            std::make_shared<ReversibleReactionConfig>(ri->id(), rj->id(), ctx)
                    );
                }
            }
        }
    }

    // look for duplicates
    for (std::size_t i = 0; i < _reversibleReactionsContainer.size(); ++i) {
        for (std::size_t j = i + 1; j < _reversibleReactionsContainer.size(); ++j) {
            const auto &revi = _reversibleReactionsContainer[i];
            const auto &revj = _reversibleReactionsContainer[j];
            if (equivalentReversibleReactions(*revi, *revj)) {
                throw std::logic_error(fmt::format(
                        "Duplicate found in reversible reactions. Registered reactions might be ambiguous "
                        "with respect to reversibility. Method: {} file: {}",
                        "DetailedBalance::searchReversibleReactions", "Actions.cpp"));
            };
        }
    }

    // create the map for lookup during runtime
    for (const auto &rev : _reversibleReactionsContainer) {
        _reversibleReactionsMap.emplace(std::make_pair(rev->forwardId, rev));
        _reversibleReactionsMap.emplace(std::make_pair(rev->backwardId, rev));
    }

    readdy::log::info(describe());
}

std::string reactions::DetailedBalance::describe() const {
    std::stringstream ss;
    if (!_reversibleReactionsContainer.empty()) {
        ss << "- reversible reactions handled subject to detailed-balance:" << readdy::util::str::newline;
        for (const auto &rev : _reversibleReactionsContainer) {
            ss << "     * -> \"" << rev->forwardReaction->name() << "\" and <- \"" << rev->backwardReaction->name()
               << "\"";
            ss << ", with types: ";
            switch (rev->reversibleType) {
                case FusionFission:
                    ss << fmt::format("{} + {} <--> {}", rev->lhsNames[0], rev->lhsNames[1], rev->rhsNames[0]);
                    break;
                case ConversionConversion:
                    ss << fmt::format("{} <--> {}", rev->lhsNames[0], rev->rhsNames[0]);
                    break;
                case EnzymaticEnzymatic:
                    ss << fmt::format("{} + {} <--> {} + {}", rev->lhsNames[0], rev->lhsNames[1],
                                      rev->rhsNames[0], rev->rhsNames[1]);
                    break;
            }
        }
    }
    return ss.str();
}

AddParticles::AddParticles(Kernel *const kernel, const std::vector<Particle> &particles)
        : particles(particles), kernel(kernel) {}

void AddParticles::perform(const util::PerformanceNode &node) {
    auto t = node.timeit();
    if(kernel != nullptr) {
        kernel->stateModel().addParticles(particles);
    } else {
        log::critical("Tried to perform {} without providing a valid kernel!", getActionName<AddParticles>());
    }
}

CalculateForces::CalculateForces() : Action() {}

top::EvaluateTopologyReactions::EvaluateTopologyReactions(scalar timeStep) : TimeStepDependentAction(timeStep) {}

}
}
}