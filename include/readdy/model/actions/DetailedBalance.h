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
 * This file contains helper classes and functions for the reaction handler 'DetailedBalance', which assures
 * detailed balance (microscopic reversibility) for reversible reactions.
 *
 * A reversible reaction `lhs <--> rhs` is defined by two reactions:
 * - forward `lhs --> rhs`
 * - backward `lhs <-- rhs`
 *
 * Supported types of reversible reactions:
 * - FusionFission A + B <--> C
 * - ConversionConversion A <--> B
 * - EnzymaticEnzymatic A + C <--> B + C, reaction radii forw. and backw. have to be equal!
 *
 * @file DetailedBalance.h
 * @brief Helpers for the reaction handler 'DetailedBalance'
 * @author chrisfroe
 * @date 05.06.18
 */

#pragma once

#include <readdy/common/common.h>
#include <readdy/model/reactions/Reaction.h>
#include <readdy/model/potentials/PotentialOrder2.h>
#include <readdy/model/Context.h>

namespace readdy {
namespace model {
namespace actions {
namespace reactions {

enum ReversibleType {
    FusionFission, ConversionConversion, EnzymaticEnzymatic
};

inline std::ostream& operator<<(std::ostream& os, const ReversibleType& reversibleType) {
    switch (reversibleType) {
        case ReversibleType::FusionFission: os << "FusionFission"; break;
        case ReversibleType::ConversionConversion: os << "ConversionConversion"; break;
        case ReversibleType::EnzymaticEnzymatic: os << "EnzymaticEnzymatic"; break;
    }
    return os;
}

struct ReversibleReactionConfig {
    using ReactionId = readdy::model::reactions::Reaction::ReactionId;
    using pot = readdy::model::potentials::PotentialOrder2 *;

    ReactionId forwardId;
    ReactionId backwardId;
    const model::reactions::Reaction *forwardReaction;
    const model::reactions::Reaction *backwardReaction;
    std::string forwardName;
    std::string backwardName;

    std::uint8_t numberLhsTypes; // either 1 or 2
    std::array<ParticleTypeId, 2> lhsTypes;
    std::array<std::string, 2> lhsNames;
    std::uint8_t numberRhsTypes; // either 1 or 2
    std::array<ParticleTypeId, 2> rhsTypes;
    std::array<std::string, 2> rhsNames;

    ReversibleType reversibleType;

    scalar microForwardRate;
    scalar microBackwardRate;

    scalar reactionRadius = 0.; // must be equal forward and backward only FusionFission and EnzymaticEnzymatic

    std::vector<pot> lhsPotentials;
    std::vector<pot> rhsPotentials; // only needed for EnzymaticEnzymatic
    scalar lhsInteractionRadius;
    scalar rhsInteractionRadius;
    scalar lhsInteractionVolume;
    scalar rhsInteractionVolume;
    scalar effectiveLhsInteractionVolume;
    scalar effectiveRhsInteractionVolume;
    scalar effectiveLhsReactionVolume;
    scalar effectiveRhsReactionVolume;

    scalar totalVolume;
    scalar kbt;
    scalar acceptancePrefactor = 1.; // EnzymaticEnzymatic only

    // To draw fission radius
    std::vector<scalar> fissionRadii; // FusionFission only
    std::vector<scalar> cumulativeFissionProb; // FusionFission only

    // these are inferred from the microscopic quantities
    scalar equilibriumConstant;
    scalar macroForwardRate;
    scalar macroBackwardRate;

    explicit ReversibleReactionConfig(ReactionId forwardId, ReactionId backwardId, const Context &ctx);

    std::string describe() const;

    scalar drawFissionDistance() const {
        auto u = readdy::model::rnd::uniform_real();
        auto it = std::lower_bound(cumulativeFissionProb.begin(), cumulativeFissionProb.end(), u);
        auto index = std::distance(cumulativeFissionProb.begin(), it);
        return fissionRadii[index];
    }
};

/**
 * Comparator for identifying equivalent reversible reactions, which are forbidden due to ambiguity.
 */
inline const bool
equivalentReversibleReactions(const ReversibleReactionConfig &rev1, const ReversibleReactionConfig &rev2) {
    return (rev1.lhsTypes == rev2.lhsTypes and rev1.rhsTypes == rev2.rhsTypes) or
           (rev1.lhsTypes == rev2.rhsTypes and rev1.rhsTypes == rev2.lhsTypes);
};

}
}
}
}
