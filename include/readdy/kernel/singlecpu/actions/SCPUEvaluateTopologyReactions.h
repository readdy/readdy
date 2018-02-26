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
 * @file SCPUEvaluateTopologyReactions.h
 * @brief << brief description >>
 * @author clonker
 * @date 18.04.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <readdy/model/actions/Actions.h>
#include <readdy/kernel/singlecpu/SCPUKernel.h>

namespace readdy {
namespace kernel {
namespace scpu {
namespace actions {
namespace top {

class SCPUEvaluateTopologyReactions : public readdy::model::actions::top::EvaluateTopologyReactions {
    using rate_t = readdy::model::top::GraphTopology::topology_reaction_rate;
public:
    SCPUEvaluateTopologyReactions(SCPUKernel* kernel, scalar timeStep);

    void perform(const util::PerformanceNode &node) override;

private:
    struct TREvent;
    using topology_reaction_events = std::vector<TREvent>;

    SCPUKernel *const kernel;

    topology_reaction_events gatherEvents();

    bool topologyDeactivated(std::size_t index) const;

    bool eventsDependent(const TREvent& evt1, const TREvent& evt2) const;

    void handleStructuralReaction(SCPUStateModel::topologies_vec &topologies,
                                  std::vector<SCPUStateModel::topology> &new_topologies,
                                  const TREvent &event, SCPUStateModel::topology_ref &topology) const;

    void handleTopologyParticleReaction(SCPUStateModel::topology_ref &topology, const TREvent &event);

    void handleTopologyTopologyReaction(SCPUStateModel::topology_ref &t1, SCPUStateModel::topology_ref &t2,
                                        const TREvent& event);

};

}
}
}
}
}
