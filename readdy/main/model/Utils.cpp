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

namespace readdy {
namespace model {
namespace util {

double getMaximumDisplacement(KernelContext& context, const double timeStep) {
    context.configure();

    double kbt = context.getKBT();

    double maximum_displacement = 0;
    for (auto &&pI : context.particleTypeRegistry().getAllRegisteredParticleTypes()) {
        double D = context.particleTypeRegistry().getDiffusionConstant(pI);
        double fMax = 0;

        for (auto &&pJ : context.particleTypeRegistry().getAllRegisteredParticleTypes()) {

            for (auto &&pot : context.potentialRegistry().getOrder2Potentials(pI, pJ)) {
                if (pot->getCutoffRadius() > 0) {
                    fMax = std::max(pot->getMaximalForce(kbt), fMax);
                }
            }
        }
        // todo magic number?
        double local_maximum_displacement = std::sqrt(2 * D * timeStep) + D * kbt * fMax * timeStep;
        maximum_displacement = std::max(local_maximum_displacement, maximum_displacement);
    }

    return maximum_displacement;
}

double getRecommendedTimeStep(unsigned int N, KernelContext& context) {
    double tau_R = 0;

    context.configure();

    double kbt = context.getKBT();
    double kReactionMax = 0;

    for (auto &&reactionO1 : context.reactionRegistry().order1_flat()) {
        kReactionMax = std::max(kReactionMax, reactionO1->getRate());
    }
    for (auto &&reactionO2 : context.reactionRegistry().order2_flat()) {
        kReactionMax = std::max(kReactionMax, reactionO2->getRate());
    }

    double tDMin = 0;
    std::unordered_map<unsigned int, double> fMaxes;
    for (auto &&pI : context.particleTypeRegistry().getAllRegisteredParticleTypes()) {
        double D = context.particleTypeRegistry().getDiffusionConstant(pI);
        double tD = 0;
        double xi = 0; // 1/(beta*Fmax)
        double fMax = 0;
        double rMin = std::numeric_limits<double>::max();

        for (auto &&reaction : context.reactionRegistry().order1_by_type(pI)) {
            if (reaction->getNProducts() == 2 && reaction->getProductDistance() > 0) {
                rMin = std::min(rMin, reaction->getProductDistance());
            }
        }

        for (auto &&pot : context.potentialRegistry().getOrder1Potentials(pI)) {
            fMax = std::max(pot->getMaximalForce(kbt), fMax);
            if (pot->getRelevantLengthScale() > 0) {
                rMin = std::min(rMin, pot->getRelevantLengthScale());
            }
        }

        for (auto &&pJ : context.particleTypeRegistry().getAllRegisteredParticleTypes()) {

            for (auto &&reaction : context.reactionRegistry().order2_by_type(pI, pJ)) {
                if (reaction->getEductDistance() > 0) {
                    rMin = std::min(rMin, reaction->getEductDistance());
                }
                if (reaction->getNProducts() == 2 && reaction->getProductDistance() > 0) {
                    rMin = std::min(rMin, reaction->getProductDistance());
                }
            }

            for (auto &&pot : context.potentialRegistry().getOrder2Potentials(pI, pJ)) {
                if (pot->getCutoffRadius() > 0) {
                    rMin = std::min(rMin, pot->getCutoffRadius());
                    fMax = std::max(pot->getMaximalForce(kbt), fMax);
                } else {
                    log::warn("Discovered potential with cutoff radius 0.");
                }
            }
        }
        double rho = rMin / 2;
        if (fMax > 0) {
            xi = 1. / (context.getKBT() * fMax);
            tD = (xi * xi / D) * (1 + rho / xi - std::sqrt(1 + 2 * rho / xi));
        } else if (D > 0) {
            tD = .5 * rho * rho / D;
        }
        fMaxes.emplace(pI, fMax);
        log::trace(" tau for {}: {} ( xi = {}, rho = {})", context.particleTypeRegistry().getParticleName(pI), tD, xi, rho);
        if (tDMin == 0) {
            tDMin = tD;
        } else {
            tDMin = std::min(tDMin, tD);
        }
    }

    log::debug("Maximal displacement for particle types per time step (stochastic + deterministic): ");
    for (auto &&pI : context.particleTypeRegistry().getAllRegisteredParticleTypes()) {
        double D = context.particleTypeRegistry().getDiffusionConstant(pI);
        double xmax = std::sqrt(2 * D * tDMin) + D * kbt * fMaxes[pI] * tDMin;
        log::debug("\t - {}: {} + {} = {}" , context.particleTypeRegistry().getParticleName(pI), std::sqrt(2 * D * tDMin),
                              D * kbt * fMaxes[pI] * tDMin, xmax);
    }

    if (kReactionMax > 0) tau_R = 1. / kReactionMax;

    double tau = std::max(tau_R, tDMin);
    if (tau_R > 0) tau = std::min(tau_R, tau);
    if (tDMin > 0) tau = std::min(tDMin, tau);
    tau /= (double) N;
    log::debug("Estimated time step: {}", tau);
    return tau;
}

}
}
}