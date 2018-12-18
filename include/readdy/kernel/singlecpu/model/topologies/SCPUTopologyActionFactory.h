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
 * @file SCPUTopologyActionFactory.h
 * @brief << brief description >>
 * @author clonker
 * @date 30.01.17
 * @copyright BSD-3
 */

#pragma once
#include <readdy/common/macros.h>
#include <readdy/model/topologies/TopologyActionFactory.h>
#include <readdy/model/topologies/reactions/TopologyReactionActionFactory.h>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(kernel)
NAMESPACE_BEGIN(scpu)
class SCPUKernel;
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(top)

namespace top = readdy::model::top;

class SCPUTopologyActionFactory : public top::TopologyActionFactory {
    SCPUKernel *const kernel;
public:
    explicit SCPUTopologyActionFactory(SCPUKernel* kernel);

    std::unique_ptr<top::pot::CalculateHarmonicBondPotential>
    createCalculateHarmonicBondPotential(const harmonic_bond * potential) const override;

    std::unique_ptr<top::pot::CalculateHarmonicAnglePotential>
    createCalculateHarmonicAnglePotential(const harmonic_angle* potential) const override;

    std::unique_ptr<top::pot::CalculateCosineDihedralPotential>
    createCalculateCosineDihedralPotential(const cos_dihedral* potential) const override;

    action_ref createChangeParticleType(top::GraphTopology* topology, const vertex &v,
                                           const ParticleTypeId &type_to) const override;

    action_ref createChangeTopologyType(top::GraphTopology * topology, const std::string &type_to) const override;

    action_ref
    createChangeParticlePosition(top::GraphTopology *topology, const vertex &v, Vec3 position) const override;

    action_ref
    createAppendParticle(top::GraphTopology *topology, const std::vector<vertex> &neighbors, ParticleTypeId type,
                         const Vec3 &position) const override;
};

NAMESPACE_END(top)
NAMESPACE_END(model)
NAMESPACE_END(scpu)
NAMESPACE_END(kernel)
NAMESPACE_END(readdy)
