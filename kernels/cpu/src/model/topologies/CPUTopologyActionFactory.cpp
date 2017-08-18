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

#include <readdy/kernel/cpu/model/topologies/CPUTopologyActionFactory.h>
#include <readdy/kernel/cpu/model/topologies/CPUTopologyActions.h>

namespace top = readdy::model::top;
namespace ctop = readdy::kernel::cpu::model::top;

ctop::CPUTopologyActionFactory::CPUTopologyActionFactory(
        readdy::kernel::cpu::CPUKernel *const kernel) : kernel(kernel) {
}

std::unique_ptr<top::pot::CalculateHarmonicBondPotential>
ctop::CPUTopologyActionFactory::createCalculateHarmonicBondPotential(
        const harmonic_bond *const potential) const {
    return std::make_unique<CPUCalculateHarmonicBondPotential>(
            &kernel->getKernelContext(), kernel->getCPUKernelStateModel().getParticleData(), potential
    );
}

std::unique_ptr<top::pot::CalculateHarmonicAnglePotential>
ctop::CPUTopologyActionFactory::createCalculateHarmonicAnglePotential(
        const harmonic_angle *const potential) const {
    return std::make_unique<CPUCalculateHarmonicAnglePotential>(
            &kernel->getKernelContext(), kernel->getCPUKernelStateModel().getParticleData(), potential
    );
}

std::unique_ptr<top::pot::CalculateCosineDihedralPotential>
ctop::CPUTopologyActionFactory::createCalculateCosineDihedralPotential(
        const cos_dihedral *const potential) const {
    return std::make_unique<CPUCalculateCosineDihedralPotential>(
            &kernel->getKernelContext(), kernel->getCPUKernelStateModel().getParticleData(), potential
    );
}

ctop::CPUTopologyActionFactory::operation_ref
ctop::CPUTopologyActionFactory::createChangeParticleType(
        top::GraphTopology *const topology, const vertex &v,
        const readdy::particle_type_type &type_to) const {
    return std::make_unique<reactions::op::CPUChangeParticleType>(kernel->getCPUKernelStateModel().getParticleData(),
                                                                  topology, v, type_to);
}
