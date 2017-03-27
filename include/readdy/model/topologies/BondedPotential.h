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
 * @file Bond.h
 * @brief << brief description >>
 * @author clonker
 * @date 26.01.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once
#include <cstddef>
#include <tuple>
#include <vector>
#include <string>
#include <readdy/model/Vec3.h>
#include "TopologyPotential.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(top)

struct BondConfiguration {
    BondConfiguration(std::size_t idx1, std::size_t idx2, double forceConstant, double length);

    std::size_t idx1, idx2;
    double length, forceConstant;
};


class BondedPotential : public TopologyPotential {
public:
    using bond_t = BondConfiguration;
    using bonds_t = std::vector<bond_t>;

    BondedPotential(Topology *const topology, const bonds_t &bonds);
    virtual ~BondedPotential() = default;

    const bonds_t &getBonds() const;
protected:
    bonds_t bonds;
};


class HarmonicBondPotential : public BondedPotential {
public:

    HarmonicBondPotential(Topology *const topology, const bonds_t &bonds);
    virtual ~HarmonicBondPotential() = default;

    double calculateEnergy(const Vec3 &x_ij, const bond_t &bond) const;

    void calculateForce(Vec3 &force, const Vec3 &x_ij, const bond_t &bond) const;

    virtual std::unique_ptr<EvaluatePotentialAction>
    createForceAndEnergyAction(const TopologyActionFactory *const) override;

};

NAMESPACE_END(top)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
