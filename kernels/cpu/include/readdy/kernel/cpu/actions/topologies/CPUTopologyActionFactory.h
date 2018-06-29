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
 * @file CPUTopologyActionFactory.h
 * @brief << brief description >>
 * @author clonker
 * @date 09.02.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once
#include <readdy/model/topologies/TopologyActionFactory.h>
#include <readdy/kernel/cpu/data/DefaultDataContainer.h>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(kernel)
NAMESPACE_BEGIN(cpu)
NAMESPACE_BEGIN(actions)
NAMESPACE_BEGIN(top)

namespace top = readdy::model::top;

class CPUTopologyActionFactory : public readdy::model::top::TopologyActionFactory {
public:
    explicit CPUTopologyActionFactory(const model::Context &context, data::DefaultDataContainer &data)
            : _context(context), _data(data) {};

    std::unique_ptr<top::pot::CalculateHarmonicBondPotential>
    createCalculateHarmonicBondPotential(const harmonic_bond* potential) const override;

    std::unique_ptr<top::pot::CalculateHarmonicAnglePotential>
    createCalculateHarmonicAnglePotential(const harmonic_angle *potential) const override;

    std::unique_ptr<top::pot::CalculateCosineDihedralPotential>
    createCalculateCosineDihedralPotential(const cos_dihedral *potential) const override;

    action_ref createChangeParticleType(top::GraphTopology* topology, const vertex &v,
                                        const ParticleTypeId &type_to) const override;

    action_ref createChangeTopologyType(top::GraphTopology *topology, const std::string &type_to) const override;

private:
    std::reference_wrapper<const model::Context> _context;
    std::reference_wrapper<data::DefaultDataContainer> _data;
};

NAMESPACE_END(top)
NAMESPACE_END(actions)
NAMESPACE_END(cpu)
NAMESPACE_END(kernel)
NAMESPACE_END(readdy)
