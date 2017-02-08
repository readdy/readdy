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
 * @file AnglePotential.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 27.01.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <readdy/model/topologies/AnglePotential.h>
#include <readdy/model/topologies/TopologyActionFactory.h>
#include <readdy/common/numeric.h>

namespace readdy {
namespace model {
namespace top {

AnglePotential::AnglePotential(Topology *const topology) : TopologyPotential(topology) {}

HarmonicAnglePotential::HarmonicAnglePotential(Topology *const topology, const angles_t &angles)
        : AnglePotential(topology), angles(angles) {}

std::unique_ptr<EvaluatePotentialAction>
HarmonicAnglePotential::createForceAndEnergyAction(const TopologyActionFactory *const factory) {
    return factory->createCalculateHarmonicAnglePotential(this);
}

const HarmonicAnglePotential::angles_t &HarmonicAnglePotential::getAngles() const {
    return angles;
}

double HarmonicAnglePotential::calculateEnergy(const Vec3 &x_ij, const Vec3 &x_kj,
                                               const HarmonicAnglePotential::Angle &angle) const {
    const double scalarProduct = x_ij * x_kj;
    const double norm_ij = std::sqrt(x_ij * x_ij);
    const double norm_kj = std::sqrt(x_kj * x_kj);
    const double theta_ijk = std::acos(scalarProduct / (norm_ij * norm_kj));
    return angle.forceConstant * (theta_ijk - angle.equilibriumAngle) * (theta_ijk - angle.equilibriumAngle);
}

void HarmonicAnglePotential::calculateForce(Vec3 &f_i, Vec3 &f_j, Vec3 &f_k, const Vec3 &x_ij, const Vec3 &x_kj,
                                            const HarmonicAnglePotential::Angle &angle) const {
    const double scalarProduct = x_ij * x_kj;
    const double norm_ij_2 = x_ij * x_ij;
    const double norm_ij = std::sqrt(norm_ij_2);
    const double norm_kj_2 = x_kj * x_kj;
    const double norm_kj = std::sqrt(norm_kj_2);
    const double norm_product = norm_ij * norm_kj;
    const double inv_norm_product = 1/norm_product;

    double cos_theta = inv_norm_product * scalarProduct;
    cos_theta = readdy::util::numeric::clamp(cos_theta, -1., 1.);

    // avoid too large values of r
    double r = std::sqrt(1.0 - cos_theta * cos_theta);
    if(r < .001) r = .001;
    r = 1./r;

    const double c = 2. * angle.forceConstant * (std::acos(cos_theta) - angle.equilibriumAngle) * r;

    const Vec3 force_i = c * cos_theta * (1/norm_ij_2) * x_ij + c * inv_norm_product * x_kj;
    const Vec3 force_k = -c * inv_norm_product * x_ij - c * cos_theta  * (1/norm_kj_2) * x_kj;

    f_i += force_i;
    f_j -= force_i + force_k;
    f_k += force_k;
}

HarmonicAnglePotential::Angle::Angle(size_t idx1, size_t idx2, size_t idx3, double theta_0, double forceConstant)
        : idx1(idx1), idx2(idx2), idx3(idx3), equilibriumAngle(theta_0), forceConstant(forceConstant) {}
}
}
}