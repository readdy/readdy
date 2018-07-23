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