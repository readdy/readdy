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
 * @file TorsionPotential.h
 * @brief << brief description >>
 * @author clonker
 * @date 26.01.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once
#include <cstddef>
#include <tuple>
#include <vector>
#include "TopologyPotential.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(top)
NAMESPACE_BEGIN(pot)

class TorsionPotential : public TopologyPotential {
public:
    TorsionPotential();
    virtual ~TorsionPotential() = default;
};

struct DihedralConfiguration {
    DihedralConfiguration(size_t idx1, size_t idx2, size_t idx3, size_t idx4, scalar forceConstant, scalar multiplicity,
             scalar equilibriumAngle);

    std::size_t idx1, idx2, idx3, idx4;
    scalar forceConstant, multiplicity, phi_0;
};

class CosineDihedralPotential : public TorsionPotential {
public:
    using dihedral_t = DihedralConfiguration;
    using dihedrals_t = std::vector<dihedral_t>;

    explicit CosineDihedralPotential(const dihedrals_t &dihedrals);
    CosineDihedralPotential(const CosineDihedralPotential&) = default;
    CosineDihedralPotential& operator=(const CosineDihedralPotential&) = default;
    CosineDihedralPotential(CosineDihedralPotential&&) = default;
    CosineDihedralPotential& operator=(CosineDihedralPotential&&) = default;
    virtual ~CosineDihedralPotential() = default;

    const dihedrals_t &getDihedrals() const;

    scalar calculateEnergy(const Vec3 &x_ji, const Vec3 &x_kj, const Vec3 &x_kl, const dihedral_t &) const;

    void calculateForce(Vec3 &f_i, Vec3 &f_j, Vec3 &f_k, Vec3 &f_l, const Vec3 &x_ji, const Vec3 &x_kj, const Vec3 &x_kl,
                        const dihedral_t &) const;

    virtual std::unique_ptr<EvaluatePotentialAction>
    createForceAndEnergyAction(const TopologyActionFactory *const factory) override;

protected:
    dihedrals_t dihedrals;
};

NAMESPACE_END(pot)
NAMESPACE_END(top)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
