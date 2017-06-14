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
 * @file SCPUEvaluateTopologyReactions.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 12.06.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <readdy/kernel/singlecpu/actions/SCPUEvaluateTopologyReactions.h>

namespace readdy {
namespace kernel {
namespace scpu {
namespace actions {
namespace top {

SCPUEvaluateTopologyReactions::SCPUEvaluateTopologyReactions(SCPUKernel *const kernel, double timeStep)
        : EvaluateTopologyReactions(timeStep), kernel(kernel) {}

void SCPUEvaluateTopologyReactions::perform() {
    using rate_t = readdy::model::top::GraphTopology::rate_t;
    auto& topologies = kernel->getSCPUKernelStateModel().topologies();

    struct TREvent {
        rate_t cumulative_rate;
        rate_t own_rate;
        std::size_t topology_idx;
        std::size_t reaction_idx;
    };

    // cumulative rate -> own rate -> topology index -> reaction index
    std::vector<TREvent> events;

    rate_t current_cumulative_rate = 0;
    std::size_t topology_idx = 0;
    for(auto& topRef : topologies) {

        if(!topRef->isDeactivated()) {

            std::size_t reaction_idx = 0;
            for (const auto &reaction : topRef->registeredReactions()) {

                TREvent event;
                event.own_rate = std::get<1>(reaction);
                event.cumulative_rate = event.own_rate + current_cumulative_rate;
                current_cumulative_rate = event.cumulative_rate;
                event.topology_idx = topology_idx;
                event.reaction_idx = reaction_idx;

                events.push_back(std::move(event));

                ++reaction_idx;
            }

        }
        ++topology_idx;
    }

}


}
}
}
}
}