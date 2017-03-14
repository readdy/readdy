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
 * @file TopologyActionFactory.h
 * @brief << brief description >>
 * @author clonker
 * @date 30.01.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once
#include <memory>
#include <readdy/common/macros.h>
#include <readdy/model/topologies/actions/TopologyActions.h>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(top)

class TopologyActionFactory {
public:
    virtual std::unique_ptr<CalculateHarmonicBondPotential>
    createCalculateHarmonicBondPotential(const HarmonicBondPotential *const) const = 0;

    virtual std::unique_ptr<CalculateHarmonicAnglePotential>
    createCalculateHarmonicAnglePotential(const HarmonicAnglePotential *const) const = 0;

    virtual std::unique_ptr<CalculateCosineDihedralPotential>
    createCalculateCosineDihedralPotential(const CosineDihedralPotential *const) const = 0;
};

NAMESPACE_END(top)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
