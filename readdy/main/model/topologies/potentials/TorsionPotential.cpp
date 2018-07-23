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
 * @file TorsionPotential.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 27.01.17
 * @copyright GPL-3
 */

#include <readdy/common/numeric.h>
#include <readdy/common/Vec3Tensor.h>
#include <readdy/model/topologies/TopologyActionFactory.h>

#define SMALL .0001

namespace readdy {
namespace model {
namespace top {
namespace pot {

scalar  CosineDihedralPotential::calculateEnergy(const Vec3 &x_ji, const Vec3 &x_kj, const Vec3 &x_kl,
                                                const dihedral_configuration &dihedral) const {
    const auto x_jk = -1. * x_kj;
    auto x_jk_norm = x_jk.norm();
    x_jk_norm = static_cast<scalar>(x_jk_norm < SMALL ? SMALL : x_jk_norm);
    const auto m = x_ji.cross(x_kj);
    const auto m_norm = m.norm();
    const auto n = x_kl.cross(x_jk);
    const auto n_norm = n.norm();
    auto m_n_norm = m_norm * n_norm;
    m_n_norm = static_cast<scalar>(m_n_norm < SMALL ? SMALL : m_n_norm);
    const scalar sin_theta = (m.cross(n) * x_jk) / (m_n_norm * x_jk_norm);
    const scalar  cos_theta = m * n / m_n_norm;
    const scalar  dih = -std::atan2(sin_theta, cos_theta);
    return dihedral.forceConstant * (1 + std::cos(dihedral.multiplicity * dih - dihedral.phi_0));
}

void
CosineDihedralPotential::calculateForce(Vec3 &f_i, Vec3 &f_j, Vec3 &f_k, Vec3 &f_l, const Vec3 &x_ji, const Vec3 &x_kj,
                                        const Vec3 &x_kl,
                                        const dihedral_configuration &dih) const {
    const auto x_jk = -1. * x_kj;
    auto x_jk_norm_squared = x_jk.normSquared();
    x_jk_norm_squared = static_cast<scalar>(x_jk_norm_squared < SMALL ? SMALL : x_jk_norm_squared);
    const auto x_jk_norm = std::sqrt(x_jk_norm_squared);
    const auto m = x_ji.cross(x_kj);
    auto m_norm_squared = m.normSquared();
    m_norm_squared = static_cast<scalar>(m_norm_squared < SMALL ? SMALL : m_norm_squared);
    const auto m_norm = std::sqrt(m_norm_squared);
    const auto n = x_kl.cross(x_jk);
    auto n_norm_squared = n.normSquared();
    n_norm_squared = static_cast<scalar>(n_norm_squared < SMALL ? SMALL : n_norm_squared);
    const auto n_norm = std::sqrt(n_norm_squared);
    auto n_m_norm = m_norm * n_norm;
    n_m_norm = static_cast<scalar>(n_m_norm < SMALL ? SMALL : n_m_norm);
    const auto m_x_n = m.cross(n);
    const auto n_m = n*m;
    const auto cos_phi = n_m / n_m_norm;
    if(std::abs(cos_phi) < SMALL) {
        return;
    }
    const auto sin_phi = (m.cross(n)) * x_jk / (n_m_norm * x_jk_norm);
    const auto phi = -std::atan2(sin_phi, cos_phi);
    const auto d_V_d_phi = -dih.forceConstant * dih.multiplicity * std::sin(dih.multiplicity * phi - dih.phi_0);
    const auto dm_norm_squared_dxi = 2 * x_jk_norm_squared * x_ji - 2 * (x_ji * x_kj) * x_kj;
    const auto dm_norm_dxi = dm_norm_squared_dxi / (2 * m_norm);
    const auto dm_norm_squared_dxk = -2 * x_ji.normSquared() * x_kj + 2 * (x_ji * x_kj) * x_ji;
    const auto dm_norm_dxk = dm_norm_squared_dxk / (2 * m_norm);
    const auto dm_norm_squared_dxj = -1. * (dm_norm_squared_dxi + dm_norm_squared_dxk);
    const auto dm_norm_dxj = dm_norm_squared_dxj / (2 * m_norm);
    // dm_squared_dxl = 0;
    const auto dn_norm_squared_dxl = 2 * x_jk_norm_squared * x_kl - 2 * (x_kl * x_jk) * x_jk;
    const auto dn_norm_dxl = dn_norm_squared_dxl / (2 * n_norm);
    const auto dn_norm_squared_dxj = -2 * x_kl.normSquared() * x_jk + 2 * (x_kl * x_jk) * x_kl;
    const auto dn_norm_dxj = dn_norm_squared_dxj / (2 * n_norm);
    const auto dn_norm_squared_dxk = -1. * (dn_norm_squared_dxl + dn_norm_squared_dxj);
    const auto dn_norm_dxk = dn_norm_squared_dxk / (2 * n_norm);

    const Vec3Tensor<3> dm_dxi{{{{0, -x_kj[2], x_kj[1]}, {x_kj[2], 0, -x_kj[0]}, {-x_kj[1], x_kj[0], 0}}}};
    const Vec3Tensor<3> dm_dxk{{{{0, -x_ji[2], x_ji[1]}, {x_ji[2], 0, -x_ji[0]}, {-x_ji[1], x_ji[0], 0}}}};
    const auto dm_dxj = (dm_dxk + dm_dxi) * -1.;

    const Vec3Tensor<3> dn_dxl = dm_dxi * -1;
    const Vec3Tensor<3> dn_dxj{{{{0, -x_kl[2], x_kl[1]}, {x_kl[2], 0, -x_kl[0]}, {-x_kl[1], x_kl[0], 0}}}};
    const auto dn_dxk = (dn_dxj + dn_dxl) * -1;

    const auto dcos_phi_prefactor = 1. / n_m_norm;
    const auto dcos_phi_dxi = dcos_phi_prefactor * (n * dm_dxi - (cos_phi * n_norm) * dm_norm_dxi);
    const auto dcos_phi_dxj = dcos_phi_prefactor * (n * dm_dxj + m * dn_dxj - (cos_phi * n_norm) * dm_norm_dxj -
                                                    (cos_phi * m_norm) * dn_norm_dxj);
    const auto dcos_phi_dxk = dcos_phi_prefactor * (n * dm_dxk + m * dn_dxk - (cos_phi * n_norm) * dm_norm_dxk -
                                                    (cos_phi * m_norm) * dn_norm_dxk);
    const auto dcos_phi_dxl = dcos_phi_prefactor * (m * dn_dxl - (cos_phi * m_norm) * dn_norm_dxl);

    const Vec3Tensor<3> dmxn_dxi{{{dm_dxi.at(0).cross(n), dm_dxi.at(1).cross(n), dm_dxi.at(2).cross(n)}}};
    const Vec3Tensor<3> dmxn_dxj{{{dm_dxj.at(0).cross(n) + m.cross(dn_dxj.at(0)),
                                          dm_dxj.at(1).cross(n) + m.cross(dn_dxj.at(1)),
                                          dm_dxj.at(2).cross(n) + m.cross(dn_dxj.at(2))}}};
    const Vec3Tensor<3> dmxn_dxk{{{dm_dxk.at(0).cross(n) + m.cross(dn_dxk.at(0)),
                                          dm_dxk.at(1).cross(n) + m.cross(dn_dxk.at(1)),
                                          dm_dxk.at(2).cross(n) + m.cross(dn_dxk.at(2))}}};
    const Vec3Tensor<3> dmxn_dxl{{{m.cross(dn_dxl.at(0)), m.cross(dn_dxl.at(1)), m.cross(dn_dxl.at(2))}}};

    const auto dsin_phi_prefactor = (1. / (n_m_norm * x_jk_norm));
    const auto dsin_phi_dxi = dsin_phi_prefactor * (x_jk * dmxn_dxi - (sin_phi * n_norm * x_jk_norm) * dm_norm_dxi);
    const auto dsin_phi_dxj =
            dsin_phi_prefactor * (x_jk * dmxn_dxj - 1. * m_x_n - (sin_phi * n_norm * x_jk_norm) * dm_norm_dxj -
                                  (sin_phi * m_norm * x_jk_norm) * dn_norm_dxj +
                                  (sin_phi * n_m_norm / x_jk_norm) * x_jk);
    const auto dsin_phi_dxk = dsin_phi_prefactor *
                              (x_jk * dmxn_dxk + m_x_n - (sin_phi * n_norm * x_jk_norm) * dm_norm_dxk -
                               (sin_phi * m_norm * x_jk_norm) * dn_norm_dxk - (sin_phi * n_m_norm / x_jk_norm) * x_jk);
    const auto dsin_phi_dxl = dsin_phi_prefactor * (x_jk * dmxn_dxl - (sin_phi * m_norm * x_jk_norm) * dn_norm_dxl);

    const auto sinphi_by_cosphi = sin_phi / cos_phi;
    const auto dphi_prefactor = -1. / (cos_phi * (1. + sinphi_by_cosphi * sinphi_by_cosphi));
    const auto dphi_dxi = dphi_prefactor * (dsin_phi_dxi - sinphi_by_cosphi * dcos_phi_dxi);
    const auto dphi_dxj = dphi_prefactor * (dsin_phi_dxj - sinphi_by_cosphi * dcos_phi_dxj);
    const auto dphi_dxk = dphi_prefactor * (dsin_phi_dxk - sinphi_by_cosphi * dcos_phi_dxk);
    const auto dphi_dxl = dphi_prefactor * (dsin_phi_dxl - sinphi_by_cosphi * dcos_phi_dxl);

    f_i += -d_V_d_phi * dphi_dxi;
    f_j += -d_V_d_phi * dphi_dxj;
    f_k += -d_V_d_phi * dphi_dxk;
    f_l += -d_V_d_phi * dphi_dxl;

}

std::unique_ptr<EvaluatePotentialAction>
CosineDihedralPotential::createForceAndEnergyAction(const TopologyActionFactory *const factory) {
    return factory->createCalculateCosineDihedralPotential(this);
}

}
}
}
}
