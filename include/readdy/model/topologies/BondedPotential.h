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

#ifndef READDY_MAIN_BOND_H
#define READDY_MAIN_BOND_H

#include <cstddef>
#include <tuple>
#include <vector>
#include <string>
#include <readdy/model/Vec3.h>
#include "TopologyPotential.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(top)

class BondedPotential : public TopologyPotential {
public:
    BondedPotential(Topology *const topology);
    virtual ~BondedPotential() = default;
};

class HarmonicBondPotential : public BondedPotential {
public:
    struct Bond;
    using bonds_t = std::vector<Bond>;

    HarmonicBondPotential(Topology *const topology, const std::vector<Bond> &bonds);
    virtual ~HarmonicBondPotential() = default;

    const bonds_t &getBonds() const;

    double calculateEnergy(const Vec3 &x_ij, const Bond &bond) const;

    void calculateForce(Vec3 &force, const Vec3 &x_ij, const Bond &bond) const;

    virtual std::unique_ptr<EvaluatePotentialAction>
    createForceAndEnergyAction(const TopologyActionFactory *const) override;

protected:
    bonds_t bonds;
};

struct HarmonicBondPotential::Bond {
    Bond(std::size_t idx1, std::size_t idx2, double forceConstant, double length);

    std::size_t idx1, idx2;
    double length, forceConstant;
};

NAMESPACE_END(top)
NAMESPACE_END(model)
NAMESPACE_END(readdy)

#endif //READDY_MAIN_BOND_H
