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
#include "TopologyPotential.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(top)

class BondPotential : public TopologyPotential {
public:
    BondPotential(Topology *const topology);
};

class HarmonicBondPotential : public BondPotential {
public:
    HarmonicBondPotential(Topology *const topology);

    void addBond(std::size_t idx1, std::size_t idx2, double length, double forceConstant);

protected:
    struct Bond;
    std::vector<Bond> bonds;
};

struct HarmonicBondPotential::Bond {
    Bond(size_t idx1, size_t idx2, double length, double forceConstant)
            : idx1(idx1), idx2(idx2), length(length), forceConstant(forceConstant) {}

    std::size_t idx1, idx2;
    double length, forceConstant;
};

NAMESPACE_END(top)
NAMESPACE_END(model)
NAMESPACE_END(readdy)

#endif //READDY_MAIN_BOND_H