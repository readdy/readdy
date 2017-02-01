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

namespace readdy {
namespace model {
namespace top {

AnglePotential::AnglePotential(Topology *const topology) : TopologyPotential(topology) {}

HarmonicAnglePotential::HarmonicAnglePotential(Topology *const topology, const angles_t &angles)
        : AnglePotential(topology), angles(angles) {}

HarmonicAnglePotential::HarmonicAnglePotential(Topology *const topology, HarmonicAnglePotential::angles_t angles)
        : AnglePotential(topology), angles(std::move(angles)) {}

std::unique_ptr<EvaluatePotentialAction>
HarmonicAnglePotential::createForceAndEnergyAction(const TopologyActionFactory *const factory) {
    return factory->createCalculateHarmonicAnglePotential(this);
}

const HarmonicAnglePotential::angles_t &HarmonicAnglePotential::getAngles() const {
    return angles;
}

double HarmonicAnglePotential::calculateForce(Vec3 &force, const Vec3 &x_ij, const Vec3 &x_kj,
                                              const HarmonicAnglePotential::Angle &angle) const {
    const double scalarProduct = x_ij * x_kj;
    const double norm_ij = std::sqrt(x_ij * x_ij);
    const double norm_kj = std::sqrt(x_kj * x_kj);
    const double theta_ijk = std::acos(scalarProduct / (norm_ij * norm_kj));
    force += (-2 * angle.forceConstant * (theta_ijk - angle.equilibriumAngle) / theta_ijk) * x_ij;
    return angle.forceConstant * (theta_ijk - angle.equilibriumAngle) * (theta_ijk - angle.equilibriumAngle);
}

double HarmonicAnglePotential::calculateEnergy(const Vec3 &x_ij, const Vec3 &x_kj,
                                               const HarmonicAnglePotential::Angle &angle) const {
    const double scalarProduct = x_ij * x_kj;
    const double norm_ij = std::sqrt(x_ij * x_ij);
    const double norm_kj = std::sqrt(x_kj * x_kj);
    const double theta_ijk = std::acos(scalarProduct / (norm_ij * norm_kj));
    return angle.forceConstant * (theta_ijk - angle.equilibriumAngle) * (theta_ijk - angle.equilibriumAngle);
}

HarmonicAnglePotential::Angle::Angle(size_t idx1, size_t idx2, size_t idx3, double theta_0, double forceConstant)
        : idx1(idx1), idx2(idx2), idx3(idx3), equilibriumAngle(theta_0), forceConstant(forceConstant) {}
}
}
}