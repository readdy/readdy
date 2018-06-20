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
 * @file Utils.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 17.11.16
 */


#include <algorithm>
#include <unordered_map>

#include <readdy/model/Utils.h>
#include <sstream>

namespace readdy {
namespace model {
namespace util {

/*scalar  getMaximumDisplacement(Context& context, const scalar  timeStep) {
    context.configure();

    scalar  kbt = context.kBT();

    scalar  maximum_displacement = 0;
    for (auto &&pI : context.particle_types().typesFlat()) {
        scalar  D = context.particle_types().diffusionConstantOf(pI);
        scalar  fMax = 0;

        for (auto &&pJ : context.particle_types().typesFlat()) {

            for (auto &&pot : context.potentials().potentialsOf(pI, pJ)) {
                if (pot->getCutoffRadius() > 0) {
                    fMax = std::max(pot->getMaximalForce(kbt), fMax);
                }
            }
        }
        // todo magic number?
        scalar  local_maximum_displacement = std::sqrt(2 * D * timeStep) + D * kbt * fMax * timeStep;
        maximum_displacement = std::max(local_maximum_displacement, maximum_displacement);
    }

    return maximum_displacement;
}

scalar  getRecommendedTimeStep(unsigned int N, Context& context) {
    scalar  tau_R = 0;

    context.configure();

    scalar  kbt = context.kBT();
    scalar  kReactionMax = 0;

    for (auto &&reactionO1 : context.reactions().order1Flat()) {
        kReactionMax = std::max(kReactionMax, reactionO1->rate());
    }
    for (auto &&reactionO2 : context.reactions().order2Flat()) {
        kReactionMax = std::max(kReactionMax, reactionO2->rate());
    }

    scalar  tDMin = 0;
    std::unordered_map<unsigned int, scalar > fMaxes;
    for (auto &&pI : context.particle_types().typesFlat()) {
        scalar  D = context.particle_types().diffusionConstantOf(pI);
        scalar  tD = 0;
        scalar  xi = 0; // 1/(beta*Fmax)
        scalar  fMax = 0;
        scalar  rMin = std::numeric_limits<scalar >::max();

        for (auto &&reaction : context.reactions().order1ByType(pI)) {
            if (reaction->nProducts() == 2 && reaction->productDistance() > 0) {
                rMin = std::min(rMin, reaction->productDistance());
            }
        }

        for (auto &&pot : context.potentials().potentialsOf(pI)) {
            fMax = std::max(pot->getMaximalForce(kbt), fMax);
            if (pot->getRelevantLengthScale() > 0) {
                rMin = std::min(rMin, pot->getRelevantLengthScale());
            }
        }

        for (auto &&pJ : context.particle_types().typesFlat()) {

            for (auto &&reaction : context.reactions().order2ByType(pI, pJ)) {
                if (reaction->eductDistance() > 0) {
                    rMin = std::min(rMin, reaction->eductDistance());
                }
                if (reaction->nProducts() == 2 && reaction->productDistance() > 0) {
                    rMin = std::min(rMin, reaction->productDistance());
                }
            }

            for (auto &&pot : context.potentials().potentialsOf(pI, pJ)) {
                if (pot->getCutoffRadius() > 0) {
                    rMin = std::min(rMin, pot->getCutoffRadius());
                    fMax = std::max(pot->getMaximalForce(kbt), fMax);
                } else {
                    log::warn("Discovered potential with cutoff radius 0.");
                }
            }
        }
        scalar  rho = rMin / 2;
        if (fMax > 0) {
            xi = static_cast<scalar>(1. / (context.kBT() * fMax));
            tD = (xi * xi / D) * (1 + rho / xi - std::sqrt(1 + 2 * rho / xi));
        } else if (D > 0) {
            tD = static_cast<scalar>(.5 * rho * rho / D);
        }
        fMaxes.emplace(pI, fMax);
        log::trace(" tau for {}: {} ( xi = {}, rho = {})", context.particle_types().nameOf(pI), tD, xi, rho);
        if (tDMin == 0) {
            tDMin = tD;
        } else {
            tDMin = std::min(tDMin, tD);
        }
    }

    log::debug("Maximal displacement for particle types per time step (stochastic + deterministic): ");
    for (auto &&pI : context.particle_types().typesFlat()) {
        scalar  D = context.particle_types().diffusionConstantOf(pI);
        scalar  xmax = std::sqrt(2 * D * tDMin) + D * kbt * fMaxes[pI] * tDMin;
        log::debug("\t - {}: {} + {} = {}" , context.particle_types().nameOf(pI), std::sqrt(2 * D * tDMin),
                              D * kbt * fMaxes[pI] * tDMin, xmax);
    }

    if (kReactionMax > 0) tau_R = static_cast<scalar>(1. / kReactionMax);

    scalar  tau = std::max(tau_R, tDMin);
    if (tau_R > 0) tau = std::min(tau_R, tau);
    if (tDMin > 0) tau = std::min(tDMin, tau);
    tau /= (scalar ) N;
    log::debug("Estimated time step: {}", tau);
    return tau;
}*/

void validateTypeName(const std::string &typeName) {
    const auto ics = invalidCharacterSequences();
    for(const auto &cs : ics) {
        if(typeName.find(cs) != std::string::npos) {
            std::stringstream ss;
            for(std::size_t i = 0; i < ics.size(); ++i) {
                if(i > 0) ss << ", ";
                ss << "\"" << ics.at(i) << "\"";
            }
            throw std::invalid_argument(fmt::format(
                    R"(Encountered invalid character sequence "{}" in type name "{}", you must not use either of {}.)",
                    cs, typeName, ss.str()
            ));
        }
    }
    if(typeName.find(' ') == 0 || typeName.find(' ') == std::string::npos-1) {
        throw std::invalid_argument(fmt::format("Type name \"{}\" contained leading/trailing whitespaces.", typeName));
    }
}

}
}
}