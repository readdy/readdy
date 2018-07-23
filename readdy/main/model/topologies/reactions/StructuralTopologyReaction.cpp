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
 * @file TopologyReaction.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 03.04.17
 * @copyright GPL-3
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
