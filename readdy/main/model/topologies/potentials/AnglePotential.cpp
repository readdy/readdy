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
 * @file AnglePotential.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 27.01.17
 * @copyright GPL-3
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
