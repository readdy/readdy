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
 * @file CPUEvaluateTopologyReactions.h
 * @brief << brief description >>
 * @author clonker
 * @date 15.06.17
 * @copyright BSD-3
 */

#pragma once


#include <readdy/model/actions/Actions.h>
#include "../CPUKernel.h"

namespace readdy::kernel::cpu::actions::top {

class CPUEvaluateTopologyReactions : public readdy::model::actions::top::EvaluateTopologyReactions {
    using rate_t = readdy::model::top::GraphTopology::ReactionRate;
public:
    CPUEvaluateTopologyReactions(CPUKernel* kernel, readdy::scalar timeStep);

    void perform() override;

private:
    struct TREvent;

    bool eventsDependent(const TREvent& evt1, const TREvent& evt2) const;

    using topology_reaction_events = std::vector<TREvent>;

    CPUKernel *const kernel;

    topology_reaction_events gatherEvents();

    bool topologyDeactivated(std::ptrdiff_t index) const;

    void handleStructuralReactionEvent(CPUStateModel::topologies_vec &topologies,
                                       std::vector<CPUStateModel::topology> &new_topologies,
                                       const TREvent &event, CPUStateModel::topology_ref &topology) const;

    void handleTopologyParticleReaction(CPUStateModel::topology_ref &topology, const TREvent &event);

    void handleTopologyTopologyReaction(CPUStateModel::topology_ref &t1, CPUStateModel::topology_ref &t2,
                                        const TREvent& event);
};


}
