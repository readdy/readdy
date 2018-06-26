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
 * @file SCPUTopologyActionFactory.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 30.01.17
 * @copyright GNU Lesser General Public License v3.0
 */


#include <readdy/kernel/singlecpu/SCPUKernel.h>

#include <readdy/kernel/singlecpu/model/topologies/SCPUTopologyActions.h>

namespace c_top = readdy::model::top;

namespace readdy {
namespace kernel {
namespace scpu {
namespace model {
namespace top {

SCPUTopologyActionFactory::SCPUTopologyActionFactory(SCPUKernel *const kernel) : kernel(kernel) {}

std::unique_ptr<c_top::pot::CalculateHarmonicBondPotential>
SCPUTopologyActionFactory::createCalculateHarmonicBondPotential(const harmonic_bond *const potential) const {
    auto& stateModel = kernel->getSCPUKernelStateModel();
    return std::make_unique<SCPUCalculateHarmonicBondPotential>(
            &kernel->context(), stateModel.getParticleData(), &stateModel.observableData(), potential
    );
}

std::unique_ptr<readdy::model::top::pot::CalculateHarmonicAnglePotential>
SCPUTopologyActionFactory::createCalculateHarmonicAnglePotential(
        const harmonic_angle *const potential) const {
    return std::make_unique<SCPUCalculateHarmonicAnglePotential>(
            &kernel->context(), kernel->getSCPUKernelStateModel().getParticleData(),
            &kernel->getSCPUKernelStateModel().observableData(), potential
    );
}

std::unique_ptr<top::pot::CalculateCosineDihedralPotential>
SCPUTopologyActionFactory::createCalculateCosineDihedralPotential(
        const cos_dihedral *const potential) const {
    return std::make_unique<SCPUCalculateCosineDihedralPotential>(
            &kernel->context(), kernel->getSCPUKernelStateModel().getParticleData(), potential
    );
}

SCPUTopologyActionFactory::action_ref
SCPUTopologyActionFactory::createChangeParticleType(top::GraphTopology *const topology, const vertex &v,
                                                    const ParticleTypeId &type_to) const {
    return std::make_unique<reactions::op::SCPUChangeParticleType>(
            kernel->getSCPUKernelStateModel().getParticleData(), topology, v, type_to
    );
}

top::reactions::actions::TopologyReactionActionFactory::action_ref
SCPUTopologyActionFactory::createChangeTopologyType(top::GraphTopology *const topology,
                                                    const std::string &type_to) const {
    return std::make_unique<readdy::model::top::reactions::actions::ChangeTopologyType>(
            topology, kernel->context().topologyRegistry().idOf(type_to)
    );
}

}
}
}
}
}
