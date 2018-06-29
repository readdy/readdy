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
 * @file TopologyReaction.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 03.04.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <readdy/model/topologies/reactions/StructuralTopologyReaction.h>

#include <readdy/model/Kernel.h>
#include <readdy/model/topologies/Utils.h>
#include <readdy/model/topologies/reactions/TopologyReactionException.h>

namespace readdy {
namespace model {
namespace top {
namespace reactions {

StructuralTopologyReaction::StructuralTopologyReaction(const reaction_function& reaction_function, const rate_function &rate_function)
        : _reaction_function(reaction_function)
        , _rate_function(rate_function) { }


StructuralTopologyReaction::StructuralTopologyReaction(const StructuralTopologyReaction::reaction_function &reaction_function, const scalar  &rate)
        : StructuralTopologyReaction(reaction_function, [rate](const GraphTopology&) -> scalar { return rate; }) {}

std::vector<GraphTopology> StructuralTopologyReaction::execute(GraphTopology &topology, const Kernel* const kernel) const {
    const auto &types = kernel->context().particleTypes();
    const auto &topology_types = kernel->context().topologyRegistry();
    auto recipe = operations(topology);
    auto& steps = recipe.steps();
    if(!steps.empty()) {
        auto topologyActionFactory = kernel->getTopologyActionFactory();
        std::vector<op::Operation::action_ptr> actions;
        actions.reserve(steps.size());
        for (auto &op : steps) {
            actions.push_back(op->create_action(&topology, topologyActionFactory));
        }

        bool exceptionOccurred = false;
        {
            // perform reaction
            auto it = actions.begin();
            try {
                for (; it != actions.end(); ++it) {
                    (*it)->execute();
                }
            } catch (const TopologyReactionException &exception) {
                log::warn("exception occurred while executing topology reaction: {}", exception.what());
                exceptionOccurred = true;
                if (raises_if_invalid()) {
                    throw;
                }
            }
            if (exceptionOccurred && rolls_back_if_invalid() && it != actions.begin()) {
                log::warn("rolling back...");
                for (;; --it) {
                    (*it)->undo();
                    if (it == actions.begin()) break;
                }
            }
        }
        if (!exceptionOccurred) {
            // post reaction
            if (expects_connected_after_reaction()) {
                bool valid = true;
                if (!topology.graph().isConnected()) {
                    // we expected it to be connected after the reaction.. but it is not, raise or rollback.
                    log::warn("The topology was expected to still be connected after the reaction, but it was not.");
                    valid = false;
                }
                {
                    // check if all particle types are topology flavored
                    for (const auto &v : topology.graph().vertices()) {
                        if (types.infoOf(v.particleType()).flavor != particleflavor::TOPOLOGY) {
                            log::warn("The topology contained particles that were not topology flavored.");
                            valid = false;
                        }
                    }
                }
                if (!valid) {
                    log::warn("GEXF representation: {}", util::to_gexf(topology.graph()));
                    if (rolls_back_if_invalid()) {
                        log::warn("rolling back...");
                        for (auto it = actions.rbegin(); it != actions.rend(); ++it) {
                            (*it)->undo();
                        }
                    } else {
                        throw TopologyReactionException(
                                "The topology was invalid after the reaction, see previous warning messages.");
                    }
                } else {
                    // if valid, update force field
                    topology.configure();
                    // and update reaction rates
                    topology.updateReactionRates(topology_types.structuralReactionsOf(topology.type()));
                }
            } else {
                if (!topology.graph().isConnected()) {
                    auto subTopologies = topology.connectedComponents();
                    assert(subTopologies.size() > 1 && "This should be at least 2 as the graph is not connected.");
                    return std::move(subTopologies);
                }
                // if valid, update force field
                topology.configure();
                // and update reaction rates
                topology.updateReactionRates(topology_types.structuralReactionsOf(topology.type()));
            }
        }
    }
    return {};
}


}
}
}
}
