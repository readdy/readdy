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
 * @file TorsionPotential.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 27.01.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <readdy/model/topologies/TorsionPotential.h>
#include <readdy/common/numeric.h>

namespace readdy {
namespace model {
namespace top {

TorsionPotential::TorsionPotential(Topology *const topology) : TopologyPotential(topology) {}

CosineDihedralPotential::CosineDihedralPotential(Topology *const t, dihedrals_t d)
        : TorsionPotential(t), dihedrals(std::move(d)) {

}

CosineDihedralPotential::CosineDihedralPotential(Topology *const topology, const dihedrals_t &dihedrals)
        : TorsionPotential(topology), dihedrals(dihedrals) {
}

const CosineDihedralPotential::dihedrals_t &CosineDihedralPotential::getDihedrals() const {
    return dihedrals;
}

double CosineDihedralPotential::calculateEnergy(const Vec3 &x_i, const Vec3 &x_j, const Vec3 &x_k, const Vec3 &x_l,
                                                const CosineDihedralPotential::Dihedral &dihedral) const {
    const auto x_ij = x_i - x_j;
    const auto x_jk = x_j - x_k;
    const auto x_lk = x_l - x_k;
    const auto m = x_ij.cross(x_jk);
    const auto m_norm = m.norm();
    const auto n = x_lk.cross(x_jk);
    const auto n_norm = n.norm();
    const double sin_theta = (n * x_ij) * x_jk.norm() / (m_norm * n_norm);
    const double cos_theta = m * n / (m_norm * n_norm);
    const double dih = -std::atan(sin_theta / cos_theta);
    return dihedral.forceConstant * (1 + std::cos(dihedral.multiplicity * dih - dihedral.equilibriumAngle));
}

CosineDihedralPotential::Dihedral::Dihedral(size_t idx1, size_t idx2, size_t idx3, size_t idx4, double forceConstant,
                                            double multiplicity, double equilibriumAngle)
        : idx1(idx1), idx2(idx2), idx3(idx3), idx4(idx4), forceConstant(forceConstant),
          equilibriumAngle(equilibriumAngle), multiplicity(multiplicity) {
    if(equilibriumAngle > readdy::util::numeric::pi() || equilibriumAngle < -readdy::util::numeric::pi()) {
        throw std::invalid_argument("the equilibrium angle should be within [-pi, pi], but was "
                                    + std::to_string(equilibriumAngle));
    }
}


}
}
}