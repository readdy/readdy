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
 * @file SCPUTopologyActionFactory.h
 * @brief << brief description >>
 * @author clonker
 * @date 30.01.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once
#include <readdy/common/macros.h>
#include <readdy/model/topologies/TopologyActionFactory.h>
#include <readdy/kernel/singlecpu/SCPUKernel.h>
#include <readdy/model/topologies/reactions/OperationFactory.h>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(kernel)
NAMESPACE_BEGIN(scpu)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(top)

namespace top = readdy::model::top;

class SCPUTopologyActionFactory : public top::TopologyActionFactory {
    const SCPUKernel *const kernel;
public:
    SCPUTopologyActionFactory(const SCPUKernel *const kernel);

    virtual std::unique_ptr<top::pot::CalculateHarmonicBondPotential>
    createCalculateHarmonicBondPotential(const harmonic_bond *const) const override;

    virtual std::unique_ptr<top::pot::CalculateHarmonicAnglePotential>
    createCalculateHarmonicAnglePotential(const harmonic_angle *const potential) const override;

    virtual std::unique_ptr<top::pot::CalculateCosineDihedralPotential>
    createCalculateCosineDihedralPotential(const cos_dihedral *const potential) const override;

    virtual operation_ref createChangeParticleType(top::GraphTopology *const topology, const vertex_t &v,
                                                   const particle_type_type &type_to) const override;
};

NAMESPACE_END(top)
NAMESPACE_END(model)
NAMESPACE_END(scpu)
NAMESPACE_END(kernel)
NAMESPACE_END(readdy)
