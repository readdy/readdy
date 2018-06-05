/********************************************************************
 * Copyright © 2018 Computational Molecular Biology Group,          *
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
    using ReactionId = readdy::model::reactions::Reaction::reaction_id;
    using pot = readdy::model::potentials::PotentialOrder2 *;

    ReactionId forwardId;
    ReactionId backwardId;
    const model::reactions::Reaction *forwardReaction;
    const model::reactions::Reaction *backwardReaction;
    std::string forwardName;
    std::string backwardName;

    std::uint8_t numberLhsTypes; // either 1 or 2
    std::array<particle_type_type, 2> lhsTypes;
    std::array<std::string, 2> lhsNames;
    std::uint8_t numberRhsTypes; // either 1 or 2
    std::array<particle_type_type, 2> rhsTypes;
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
