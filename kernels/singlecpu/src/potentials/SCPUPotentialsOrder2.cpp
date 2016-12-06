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
 * SingleCPU implementations of potentials of second order (particle-particle interactions).
 * This is where actual forces and energies are calculated.
 *
 * @file PotentialsOrder2.cpp
 * @brief Calculate forces and energies of second-order potentials.
 * @author clonker
 * @date 09.06.16
 */

#include <readdy/kernel/singlecpu/potentials/SCPUPotentialsOrder2.h>
#include <readdy/kernel/singlecpu/SCPUKernel.h>

namespace readdy {
namespace kernel {
namespace scpu {
namespace potentials {

double SCPUHarmonicRepulsion::calculateEnergy(const vec_t &x_ij) const {
    auto distanceSquared = x_ij * x_ij;
    if (distanceSquared < getSumOfParticleRadiiSquared()) {
        distanceSquared = std::sqrt(distanceSquared);
        distanceSquared -= getSumOfParticleRadii();
        distanceSquared *= distanceSquared;
        return 0.5 * distanceSquared * getForceConstant();
    } else {
        return 0;
    }
}

void SCPUHarmonicRepulsion::calculateForce(vec_t &force, const vec_t &x_ij) const {
    auto squared = x_ij * x_ij;
    if (squared < getSumOfParticleRadiiSquared() && squared > 0) {
        squared = std::sqrt(squared);
        force = (getForceConstant() * (squared - getSumOfParticleRadii())) / squared * x_ij;
    } else {
        force = {0, 0, 0};
    }
}

void SCPUHarmonicRepulsion::calculateForceAndEnergy(vec_t &force, double &energy, const vec_t &x_ij) const {
    auto squared = x_ij * x_ij;
    if (squared < getSumOfParticleRadiiSquared() && squared > 0) {
        squared = std::sqrt(squared);
        energy += 0.5 * getForceConstant() * std::pow(squared - getSumOfParticleRadii(), 2);
        force = (getForceConstant() * (squared - getSumOfParticleRadii())) / squared * x_ij;
    } else {
        force = {0, 0, 0};
    }
}

SCPUHarmonicRepulsion::SCPUHarmonicRepulsion(const readdy::model::Kernel *const kernel)
        : readdy::model::potentials::HarmonicRepulsion(kernel) {}

double SCPUHarmonicRepulsion::getCutoffRadius() const {
    return sumOfParticleRadii;
}

double SCPUHarmonicRepulsion::getCutoffRadiusSquared() const {
    return sumOfParticleRadiiSquared;
}


SCPUWeakInteractionPiecewiseHarmonic::SCPUWeakInteractionPiecewiseHarmonic(const readdy::model::Kernel *const kernel)
        : readdy::model::potentials::WeakInteractionPiecewiseHarmonic(kernel) {}

double SCPUWeakInteractionPiecewiseHarmonic::calculateEnergy(const readdy::model::Vec3 &x_ij) const {
    const auto dist = std::sqrt(x_ij * x_ij);
    const auto len_part2 = noInteractionDistance - desiredParticleDistance;
    if (dist < desiredParticleDistance) {
        // repulsive as we are closer than the desired distance
        return .5 * forceConstant * (dist - desiredParticleDistance) * (dist - desiredParticleDistance) -
               depthAtDesiredDistance;
    } else {
        // attractive as we are further (but not too far) apart than the desired distance
        if (dist < desiredParticleDistance + .5 * len_part2) {
            return .5 * depthAtDesiredDistance * (1 / (.5 * len_part2)) * (1 / (.5 * len_part2)) *
                   (dist - desiredParticleDistance) * (dist - desiredParticleDistance) -
                   depthAtDesiredDistance;
        } else {
            // if we are not too far apart but still further than in the previous case, attractive
            if (dist < noInteractionDistance) {
                return -0.5 * depthAtDesiredDistance * (1 / (0.5 * len_part2)) * (1 / (0.5 * len_part2)) *
                       (dist - noInteractionDistance) * (dist - noInteractionDistance);
            }
        }
    }
    return 0;
}

void
SCPUWeakInteractionPiecewiseHarmonic::calculateForce(readdy::model::Vec3 &force, const readdy::model::Vec3 &x_ij) const {
    const auto dist = std::sqrt(x_ij * x_ij);
    const auto len_part2 = noInteractionDistance - desiredParticleDistance;
    double factor = 0;
    if (dist < desiredParticleDistance) {
        // repulsive as we are closer than the desired distance
        factor = -1 * forceConstant * (desiredParticleDistance - dist);
    } else {
        // attractive as we are further (but not too far) apart than the desired distance
        if (dist < desiredParticleDistance + .5 * len_part2) {
            factor = -1 * depthAtDesiredDistance * (1 / (.5 * len_part2)) * (1 / (.5 * len_part2)) *
                     (desiredParticleDistance - dist);
        } else {
            // if we are not too far apart but still further than in the previous case, attractive
            if (dist < noInteractionDistance) {
                factor = depthAtDesiredDistance * (1 / (.5 * len_part2)) * (1 / (.5 * len_part2)) *
                         (noInteractionDistance - dist);
            }
        }
    }
    if (dist > 0 && factor != 0) {
        force = factor * x_ij / dist;
    } else {
        force = {0, 0, 0};
    }
}

void SCPUWeakInteractionPiecewiseHarmonic::calculateForceAndEnergy(readdy::model::Vec3 &force, double &energy,
                                                               const readdy::model::Vec3 &x_ij) const {
    energy += calculateEnergy(x_ij);
    calculateForce(force, x_ij);
}

double SCPUWeakInteractionPiecewiseHarmonic::getCutoffRadius() const {
    return noInteractionDistance;
}

double SCPUWeakInteractionPiecewiseHarmonic::getCutoffRadiusSquared() const {
    return noInteractionDistanceSquared;
}


}
}
}
}

