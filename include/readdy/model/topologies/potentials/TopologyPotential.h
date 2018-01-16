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
 * @file TopologyPotential.h
 * @brief << brief description >>
 * @author clonker
 * @date 27.01.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once
#include <readdy/common/macros.h>
#include <readdy/model/topologies/potentials/TopologyPotentialAction.h>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(top)
class Topology;
class TopologyActionFactory;
NAMESPACE_BEGIN(pot)

class TopologyPotential {
public:
    TopologyPotential() = default;

    TopologyPotential(const TopologyPotential&) = default;
    TopologyPotential& operator=(const TopologyPotential&) = delete;
    TopologyPotential(TopologyPotential&&) = default;
    TopologyPotential& operator=(TopologyPotential&&) = delete;

    virtual ~TopologyPotential() = default;

    virtual std::unique_ptr<EvaluatePotentialAction> createForceAndEnergyAction(const TopologyActionFactory*) = 0;
};

NAMESPACE_END(pot)
NAMESPACE_END(top)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
