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
 * @file AnglePotential.h
 * @brief << brief description >>
 * @author clonker
 * @date 26.01.17
 * @copyright GNU Lesser General Public License v3.0
 */

#ifndef READDY_MAIN_ANGLEPOTENTIAL_H
#define READDY_MAIN_ANGLEPOTENTIAL_H

#include <cstddef>
#include <tuple>
#include <vector>
#include "TopologyPotential.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(top)

class AnglePotential : public TopologyPotential{
public:
    using angles_t = std::vector<std::tuple<std::size_t, std::size_t, std::size_t>>;
    AnglePotential(Topology *const topology);
};

class HarmonicAnglePotential : public AnglePotential{
public:
    struct Angle;
    using angles_t = std::vector<Angle>;
    HarmonicAnglePotential(Topology *const topology, const angles_t& angles);
    HarmonicAnglePotential(Topology *const topology, angles_t angles);

    virtual std::unique_ptr<EvaluatePotentialAction>
    createForceAndEnergyAction(const TopologyActionFactory *const factory) override;

    const angles_t &getAngles() const;

    double calculateEnergy(const Vec3& x_ij, const Vec3& x_kj, const Angle& angle) const;

    double calculateForce(Vec3 &force, const Vec3 &x_ij, const Vec3 &x_kj, const Angle &angle) const;

protected:
    angles_t angles;
};

struct HarmonicAnglePotential::Angle {

    Angle(size_t idx1, size_t idx2, size_t idx3, double equilibriumAngle, double forceConstant);

    const std::size_t idx1, idx2, idx3;
    const double equilibriumAngle, forceConstant;
};

NAMESPACE_END(top)
NAMESPACE_END(model)
NAMESPACE_END(readdy)

#endif //READDY_MAIN_ANGLEPOTENTIAL_H
