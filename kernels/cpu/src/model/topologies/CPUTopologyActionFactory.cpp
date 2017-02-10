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

readdy::kernel::cpu::model::top::CPUTopologyActionFactory::CPUTopologyActionFactory(
        readdy::kernel::cpu::CPUKernel *const kernel) : kernel(kernel) {
}

std::unique_ptr<readdy::model::top::CalculateHarmonicBondPotential>
readdy::kernel::cpu::model::top::CPUTopologyActionFactory::createCalculateHarmonicBondPotential(
        const readdy::model::top::HarmonicBondPotential *const potential) const {
    return std::make_unique<CPUCalculateHarmonicBondPotential>(
            &kernel->getKernelContext(), kernel->getCPUKernelStateModel().getParticleData(), potential
    );
}

std::unique_ptr<readdy::model::top::CalculateHarmonicAnglePotential>
readdy::kernel::cpu::model::top::CPUTopologyActionFactory::createCalculateHarmonicAnglePotential(
        const readdy::model::top::HarmonicAnglePotential *const potential) const {
    return std::make_unique<CPUCalculateHarmonicAnglePotential>(
            &kernel->getKernelContext(), kernel->getCPUKernelStateModel().getParticleData(), potential
    );
}

std::unique_ptr<readdy::model::top::CalculateCosineDihedralPotential>
readdy::kernel::cpu::model::top::CPUTopologyActionFactory::createCalculateCosineDihedralPotential(
        const readdy::model::top::CosineDihedralPotential *const potential) const {
    return std::make_unique<CPUCalculateCosineDihedralPotential>(
            &kernel->getKernelContext(), kernel->getCPUKernelStateModel().getParticleData(), potential
    );
}
