/********************************************************************
 * Copyright © 2018 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * Redistribution and use in source and binary forms, with or       *
 * without modification, are permitted provided that the            *
 * following conditions are met:                                    *
 *  1. Redistributions of source code must retain the above         *
 *     copyright notice, this list of conditions and the            *
 *     following disclaimer.                                        *
 *  2. Redistributions in binary form must reproduce the above      *
 *     copyright notice, this list of conditions and the following  *
 *     disclaimer in the documentation and/or other materials       *
 *     provided with the distribution.                              *
 *  3. Neither the name of the copyright holder nor the names of    *
 *     its contributors may be used to endorse or promote products  *
 *     derived from this software without specific                  *
 *     prior written permission.                                    *
 *                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND           *
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,      *
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF         *
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE         *
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR            *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,     *
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,         *
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER *
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,      *
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)    *
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF      *
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
 ********************************************************************/


/**
 * << detailed description >>
 *
 * @file SCPUTopologyActionFactory.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 30.01.17
 * @copyright GPL-3
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

top::reactions::actions::TopologyReactionActionFactory::action_ref
SCPUTopologyActionFactory::createChangeParticlePosition(top::GraphTopology *topology,
                                                        const vertex &v,
                                                        Vec3 position) const {
    return std::make_unique<reactions::op::SCPUChangeParticlePosition>(
            kernel->getSCPUKernelStateModel().getParticleData(), topology, v, position
    );
}

}
}
}
}
}
