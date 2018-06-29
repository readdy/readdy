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
 * @file CPUTopologyActionFactory.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 09.02.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <readdy/kernel/cpu/actions/topologies/CPUTopologyActionFactory.h>
#include <readdy/kernel/cpu/actions/topologies/CPUTopologyActions.h>

namespace top = readdy::model::top;
namespace ctop = readdy::kernel::cpu::actions::top;

std::unique_ptr<top::pot::CalculateHarmonicBondPotential>
ctop::CPUTopologyActionFactory::createCalculateHarmonicBondPotential(
        const harmonic_bond *const potential) const {
    return std::make_unique<CPUCalculateHarmonicBondPotential>(
            &_context.get(), &_data.get(), potential
    );
}

std::unique_ptr<top::pot::CalculateHarmonicAnglePotential>
ctop::CPUTopologyActionFactory::createCalculateHarmonicAnglePotential(
        const harmonic_angle *const potential) const {
    return std::make_unique<CPUCalculateHarmonicAnglePotential>(
            &_context.get(), &_data.get(), potential
    );
}

std::unique_ptr<top::pot::CalculateCosineDihedralPotential>
ctop::CPUTopologyActionFactory::createCalculateCosineDihedralPotential(
        const cos_dihedral *const potential) const {
    return std::make_unique<CPUCalculateCosineDihedralPotential>(
            &_context.get(), &_data.get(), potential
    );
}

ctop::CPUTopologyActionFactory::action_ref
ctop::CPUTopologyActionFactory::createChangeParticleType(
        top::GraphTopology *const topology, const vertex &v,
        const readdy::ParticleTypeId &type_to) const {
    return std::make_unique<reactions::op::CPUChangeParticleType>(&_data.get(), topology, v, type_to);
}

top::reactions::actions::TopologyReactionActionFactory::action_ref
readdy::kernel::cpu::actions::top::CPUTopologyActionFactory::createChangeTopologyType(top::GraphTopology *const topology,
                                                                                    const std::string &type_to) const {
    return std::make_unique<readdy::model::top::reactions::actions::ChangeTopologyType>(
            topology, _context.get().topologyRegistry().idOf(type_to)
    );
}
