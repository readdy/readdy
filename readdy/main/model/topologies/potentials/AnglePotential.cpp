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

#include <readdy/common/numeric.h>
#include <readdy/model/topologies/TopologyActionFactory.h>

#define SMALL 1e-10

namespace readdy {
namespace model {
namespace top {
namespace pot {

std::unique_ptr<EvaluatePotentialAction>
HarmonicAnglePotential::createForceAndEnergyAction(const TopologyActionFactory *const factory) {
    return factory->createCalculateHarmonicAnglePotential(this);
}

scalar HarmonicAnglePotential::calculateEnergy(const Vec3 &x_ij, const Vec3 &x_kj,
                                               const angle &angle) const {
    const scalar scalarProduct = x_ij * x_kj;
    const scalar norm_ij = x_ij.norm();
    const scalar norm_kj = x_kj.norm();
    const scalar theta_ijk = std::acos(scalarProduct / (norm_ij * norm_kj));
    return angle.forceConstant * (theta_ijk - angle.equilibriumAngle) * (theta_ijk - angle.equilibriumAngle);
}

void HarmonicAnglePotential::calculateForce(Vec3 &f_i, Vec3 &f_j, Vec3 &f_k, const Vec3 &x_ji, const Vec3 &x_jk,
                                            const angle &angle) const {
    const scalar scalarProduct = x_ji * x_jk;
    scalar norm_ji_2 = x_ji * x_ji;
    if (norm_ji_2 < SMALL) {
        norm_ji_2 = SMALL;
    }
    const scalar norm_ji = std::sqrt(norm_ji_2);
    scalar norm_jk_2 = x_jk * x_jk;
    if (norm_jk_2 < SMALL) {
        norm_jk_2 = SMALL;
    }
    const scalar norm_jk = std::sqrt(norm_jk_2);
    scalar norm_product = norm_ji * norm_jk;
    if (norm_product < SMALL) {
        norm_product = SMALL;
    }
    const scalar inv_norm_product = c_::one / norm_product;

    scalar cos_theta = inv_norm_product * scalarProduct;
    cos_theta = readdy::util::numeric::clamp(cos_theta, -c_::one, c_::one);

    scalar sin_theta_inv = std::sqrt(c_::one - cos_theta * cos_theta);
    // avoid too small values of sin_theta
    if (sin_theta_inv < SMALL) {
        sin_theta_inv = SMALL;
    }
    sin_theta_inv = c_::one / sin_theta_inv;

    const scalar c = c_::two * angle.forceConstant * (std::acos(cos_theta) - angle.equilibriumAngle) * sin_theta_inv;

    const Vec3 force_i = c * cos_theta * (c_::one / norm_ji_2) * x_ji - c * inv_norm_product * x_jk;
    const Vec3 force_k = -c * inv_norm_product * x_ji + c * cos_theta * (c_::one / norm_jk_2) * x_jk;

    f_i -= force_i;
    f_j += force_i + force_k;
    f_k -= force_k;
}

}
}
}
}
