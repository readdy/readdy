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
 * @file TopologyActions.h
 * @brief << brief description >>
 * @author clonker
 * @date 30.01.17
 * @copyright GNU Lesser General Public License v3.0
 */

#ifndef READDY_MAIN_TOPOLOGYACTIONS_H
#define READDY_MAIN_TOPOLOGYACTIONS_H

#include "TopologyAction.h"
#include <readdy/model/topologies/BondedPotential.h>
#include <readdy/model/topologies/AnglePotential.h>
#include <readdy/model/topologies/TorsionPotential.h>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(top)

class CalculateHarmonicBondPotential : public EvaluatePotentialAction {
public:
    CalculateHarmonicBondPotential(const KernelContext *const context) : EvaluatePotentialAction(context) {}
};

class CalculateHarmonicAnglePotential : public EvaluatePotentialAction {
public:
    CalculateHarmonicAnglePotential(const KernelContext *const context) : EvaluatePotentialAction(context) {}
};

class CalculateCosineDihedralPotential : public EvaluatePotentialAction {
public:
    CalculateCosineDihedralPotential(const KernelContext *const context) : EvaluatePotentialAction(context) {}
};

NAMESPACE_END(top)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
#endif //READDY_MAIN_TOPOLOGYACTIONS_H
