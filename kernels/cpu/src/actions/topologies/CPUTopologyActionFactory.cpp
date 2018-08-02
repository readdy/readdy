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
 * @file CPUTopologyActionFactory.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 09.02.17
 * @copyright GPL-3
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

top::reactions::actions::TopologyReactionActionFactory::action_ref
readdy::kernel::cpu::actions::top::CPUTopologyActionFactory::createChangeParticlePosition(
        top::GraphTopology *topology, const top::reactions::actions::TopologyReactionActionFactory::vertex &v,
        readdy::Vec3 position) const {
    return std::make_unique<reactions::op::CPUChangeParticlePosition>(&_data.get(), topology, v, position);
}
