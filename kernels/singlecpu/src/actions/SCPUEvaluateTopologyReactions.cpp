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

/**
 * Struct holding information about a topology reaction event.
 * If:
 *   - reaction_idx == -1: external topology reaction
 *   - reaction_idx >= 0: internal (fission/conversion-type) topology reaction
 */
struct SCPUEvaluateTopologyReactions::TREvent {
    using index_type = std::size_t;

    rate_t cumulative_rate{0};
    rate_t own_rate{0};
    std::size_t topology_idx{0};
    std::ptrdiff_t topology_idx2{-1};
    std::size_t reaction_idx{0};
    particle_type_type t1{0}, t2{0};
    // idx1 is always the particle that belongs to a topology
    index_type idx1{0}, idx2{0};
    bool spatial {false};

};


SCPUEvaluateTopologyReactions::SCPUEvaluateTopologyReactions(SCPUKernel *const kernel, scalar timeStep)
        : EvaluateTopologyReactions(timeStep), kernel(kernel) {}

template<bool approximated>
bool performReactionEvent(const scalar rate, const scalar timeStep) {
    if (approximated) {
        return readdy::model::rnd::uniform_real() < rate * timeStep;
    } else {
        return readdy::model::rnd::uniform_real() < 1 - std::exp(-rate * timeStep);
    }
}

bool shouldPerformEvent(const scalar rate, const scalar timeStep, bool approximated) {
    return approximated ? performReactionEvent<true>(rate, timeStep) : performReactionEvent<false>(rate, timeStep);
}

void SCPUEvaluateTopologyReactions::perform(const util::PerformanceNode &node) {
    auto t = node.timeit();
    auto &model = kernel->getSCPUKernelStateModel();
    const auto &context = kernel->getKernelContext();
    auto &topologies = model.topologies();

    if (!topologies.empty()) {

        auto events = gatherEvents();

        if (!events.empty()) {
            auto end = events.end();

            std::vector<readdy::model::top::GraphTopology> new_topologies;
            while (end != events.begin()) {
                const auto cumulative_rate = (end - 1)->cumulative_rate;

                const auto x = readdy::model::rnd::uniform_real(c_::zero, cumulative_rate);

                const auto eventIt = std::lower_bound(
                        events.begin(), end, x, [](const TREvent &elem1, const rate_t elem2) {
                            return elem1.cumulative_rate < elem2;
                        }
                );

                if (eventIt != events.end()) {
                    const auto &event = *eventIt;

                    if (performReactionEvent<true>(event.own_rate, timeStep)) {
                        log::trace("picked event {} / {} with rate {}", std::distance(events.begin(), eventIt) + 1,
                                   events.size(), eventIt->own_rate);
                        // perform the event!
                        auto &topology = topologies.at(event.topology_idx);
                        if (topology->isDeactivated()) {
                            log::critical("deactivated topology with idx {}", event.topology_idx);
                            for (auto it = events.begin(); it != end; ++it) {
                                log::warn(" -> event {}, {}, {}, {}", it->topology_idx, it->reaction_idx, it->own_rate,
                                          it->cumulative_rate);
                            }
                            for (auto it = end; it != events.end(); ++it) {
                                log::warn(" -> deactivated event {}, {}, {}, {}", it->topology_idx, it->reaction_idx,
                                          it->own_rate, it->cumulative_rate);
                            }
                        }
                        assert(!topology->isDeactivated());
                        if(!event.spatial) {
                            handleStructuralReaction(topologies, new_topologies, event, topology);
                        } else {
                            if(event.topology_idx2 >= 0) {
                                auto &top2 = topologies.at(static_cast<std::size_t>(event.topology_idx2));
                                handleTopologyTopologyReaction(topology, top2, event);
                            } else {
                                handleTopologyParticleReaction(topology, event);
                            }
                        }

                        // loop over all other events, swap events considering this particular topology to the end
                        rate_t cumsum = 0;
                        //log::warn("------");
                        for (auto it = events.begin();;) {
                            if (eventsDependent(*it, event)) {
                                --end;
                                //log::warn("swapping event with topology idx {}", it->topology_idx);
                                std::iter_swap(it, end);
                            }
                            cumsum += it->own_rate;
                            it->cumulative_rate = cumsum;
                            if (it >= end) {
                                break;
                            }
                            if (!eventsDependent(*it, event)) {
                                ++it;
                            }
                        }

                    } else {

                        // remove event from the list (ie shift it to the end)
                        --end;

                        // swap with last element that is not yet deactivated
                        std::iter_swap(eventIt, end);

                        // this element's cumulative rate gets initialized with its own rate
                        eventIt->cumulative_rate = eventIt->own_rate;
                        if (eventIt > events.begin()) {
                            // and then increased by the cumulative rate of its predecessor
                            eventIt->cumulative_rate += (eventIt - 1)->cumulative_rate;
                        }
                        // now update the cumulative rates of all following elements...
                        auto cumsum = (*eventIt).cumulative_rate;
                        for (auto _it = eventIt + 1; _it < end; ++_it) {
                            cumsum += (*_it).own_rate;
                            (*_it).cumulative_rate = cumsum;
                        }

                    }

                } else {
                    log::critical("this should not happen (event not found by drawn rate)");
                    throw std::logic_error("this should not happen (event not found by drawn rate)");
                }
            }

            if (!new_topologies.empty()) {
                for (auto &&top : new_topologies) {
                    if (!top.isNormalParticle(*kernel)) {
                        // we have a new topology here, update data accordingly.
                        top.updateReactionRates(context.topology_registry().structural_reactions_of(top.type()));
                        top.configure();
                        model.insert_topology(std::move(top));
                    } else {
                        // if we have a single particle that is not of flavor topology, remove from topology structure!
                        model.getParticleData()->entry_at(top.getParticles().front()).topology_index = -1;
                    }
                }
            }
        }
    }

}

SCPUEvaluateTopologyReactions::topology_reaction_events SCPUEvaluateTopologyReactions::gatherEvents() {
    const auto& context = kernel->getKernelContext();
    const auto& topology_registry = context.topology_registry();
    topology_reaction_events events;
    {
        rate_t current_cumulative_rate = 0;
        std::size_t topology_idx = 0;
        for (auto &top : kernel->getSCPUKernelStateModel().topologies()) {
            if (!top->isDeactivated()) {
                std::size_t reaction_idx = 0;
                for (const auto &reaction : topology_registry.structural_reactions_of(top->type())) {
                    TREvent event{};
                    event.own_rate = top->rates().at(reaction_idx);
                    event.cumulative_rate = event.own_rate + current_cumulative_rate;
                    current_cumulative_rate = event.cumulative_rate;
                    event.topology_idx = topology_idx;
                    event.reaction_idx = reaction_idx;

                    events.push_back(event);
                    ++reaction_idx;
                }
            }
            ++topology_idx;
        }


        if (!topology_registry.spatial_reaction_registry().empty()) {
            const auto &d2 = context.getDistSquaredFun();
            const auto &stateModel = kernel->getSCPUKernelStateModel();
            const auto &data = *stateModel.getParticleData();
            const auto &nl = *stateModel.getNeighborList();

            for (auto pair : nl) {
                auto &entry = data.entry_at(pair.idx1);
                auto &neighbor = data.entry_at(pair.idx2);
                if (!entry.deactivated && !neighbor.deactivated) {
                    topology_type_type tt1 = entry.topology_index >= 0 ? stateModel.topologies().at(static_cast<std::size_t>(entry.topology_index))->type() : static_cast<topology_type_type>(-1);
                    topology_type_type tt2 = neighbor.topology_index >= 0 ? stateModel.topologies().at(static_cast<std::size_t>(neighbor.topology_index))->type() : static_cast<topology_type_type>(-1);
                    if(tt1 == -1 && tt2 == -1) continue;
                    const auto &reactions = topology_registry.spatial_reactions_by_type(entry.type, tt1,
                                                                                        neighbor.type, tt2);
                    const auto distSquared = d2(entry.pos, neighbor.pos);
                    std::size_t reaction_index = 0;
                    for (const auto &reaction : reactions) {
                        if (distSquared < reaction.radius() * reaction.radius()) {
                            TREvent event{};
                            event.own_rate = reaction.rate();
                            event.cumulative_rate = event.own_rate + current_cumulative_rate;
                            current_cumulative_rate = event.cumulative_rate;
                            switch (reaction.mode()) {
                                case readdy::model::top::reactions::STRMode::TT_FUSION:
                                    if(entry.topology_index >= 0 && neighbor.topology_index >= 0 && entry.topology_index != neighbor.topology_index) {
                                        event.topology_idx = static_cast<std::size_t>(entry.topology_index);
                                        event.topology_idx2 = neighbor.topology_index;
                                        event.t1 = entry.type;
                                        event.t2 = neighbor.type;
                                        event.idx1 = pair.idx1;
                                        event.idx2 = pair.idx2;
                                        event.reaction_idx = reaction_index;
                                        event.spatial = true;

                                        events.push_back(event);
                                    }
                                    break;
                                case readdy::model::top::reactions::STRMode::TT_ENZYMATIC: // fall through
                                case readdy::model::top::reactions::STRMode::TT_FUSION_ALLOW_SELF:
                                    if(entry.topology_index >= 0 && neighbor.topology_index >= 0) {
                                        event.topology_idx = static_cast<std::size_t>(entry.topology_index);
                                        event.topology_idx2 = neighbor.topology_index;
                                        event.t1 = entry.type;
                                        event.t2 = neighbor.type;
                                        event.idx1 = pair.idx1;
                                        event.idx2 = pair.idx2;
                                        event.reaction_idx = reaction_index;
                                        event.spatial = true;

                                        events.push_back(event);
                                    }
                                    break;
                                case readdy::model::top::reactions::STRMode::TP_ENZYMATIC: // fall through
                                case readdy::model::top::reactions::STRMode::TP_FUSION:
                                    if (entry.topology_index >= 0 && neighbor.topology_index < 0) {
                                        event.topology_idx = static_cast<std::size_t>(entry.topology_index);
                                        event.t1 = entry.type;
                                        event.t2 = neighbor.type;
                                        event.idx1 = pair.idx1;
                                        event.idx2 = pair.idx2;
                                        event.reaction_idx = reaction_index;
                                        event.spatial = true;

                                        events.push_back(event);
                                    } else if (entry.topology_index < 0 && neighbor.topology_index >= 0) {
                                        event.topology_idx = static_cast<std::size_t>(neighbor.topology_index);
                                        event.t1 = neighbor.type;
                                        event.t2 = entry.type;
                                        event.idx1 = pair.idx2;
                                        event.idx2 = pair.idx1;
                                        event.reaction_idx = reaction_index;
                                        event.spatial = true;

                                        events.push_back(event);
                                    }
                                    break;
                            }
                        }
                        ++reaction_index;
                    }
                }
            }
        }
    }
    return events;
}

bool SCPUEvaluateTopologyReactions::eventsDependent(const SCPUEvaluateTopologyReactions::TREvent &evt1,
                                                    const SCPUEvaluateTopologyReactions::TREvent &evt2) const {
    if(evt1.topology_idx == evt2.topology_idx || evt1.topology_idx == evt2.topology_idx2) {
        return true;
    }
    return evt1.topology_idx2 >= 0 && (evt1.topology_idx2 == evt2.topology_idx || evt1.topology_idx2 == evt2.topology_idx2);
}

void SCPUEvaluateTopologyReactions::handleTopologyTopologyReaction(SCPUStateModel::topology_ref &t1,
                                                                   SCPUStateModel::topology_ref &t2,
                                                                   const SCPUEvaluateTopologyReactions::TREvent &event) {
    const auto& context = kernel->getKernelContext();
    const auto& top_registry = context.topology_registry();
    const auto& reaction = top_registry.spatial_reactions_by_type(event.t1, t1->type(), event.t2, t2->type()).at(event.reaction_idx);

    auto& model = kernel->getSCPUKernelStateModel();
    auto& data = *model.getParticleData();

    auto& entry1 = data.entry_at(event.idx1);
    auto& entry2 = data.entry_at(event.idx2);
    auto& entry1Type = entry1.type;
    auto& entry2Type = entry2.type;
    auto top_type_to1 = reaction.top_type_to1();
    auto top_type_to2 = reaction.top_type_to2();
    if(entry1Type == reaction.type1() && t1->type() == reaction.top_type1()) {
        entry1Type = reaction.type_to1();
        entry2Type = reaction.type_to2();
    } else {
        std::swap(top_type_to1, top_type_to2);
        entry1Type = reaction.type_to2();
        entry2Type = reaction.type_to1();
    }

    // topology - topology fusion
    if(reaction.is_fusion()) {
        //topology->appendTopology(otherTopology, ...)
        if(event.topology_idx == event.topology_idx2) {
            // introduce edge if not already present
            auto v1 = t1->vertexForParticle(event.idx1);
            auto v2 = t1->vertexForParticle(event.idx2);
            if(!t1->graph().containsEdge(v1, v2)) {
                t1->graph().addEdge(v1, v2);
            }
        } else {
            // merge topologies
            for(auto pidx : t2->getParticles()) {
                data.entry_at(pidx).topology_index = event.topology_idx;
            }
            auto &topologies = model.topologies();
            t1->appendTopology(*t2, event.idx2, entry2Type, event.idx1, entry1Type, top_type_to1);
            topologies.erase(topologies.begin() + event.topology_idx2);
        }
    } else {
        t1->vertexForParticle(event.idx1)->setParticleType(entry1Type);
        t2->vertexForParticle(event.idx2)->setParticleType(entry2Type);
        t1->type() = top_type_to1;
        t2->type() = top_type_to2;

        t2->updateReactionRates(context.topology_registry().structural_reactions_of(t2->type()));
        t2->configure();
    }
    t1->updateReactionRates(context.topology_registry().structural_reactions_of(t1->type()));
    t1->configure();
}

void SCPUEvaluateTopologyReactions::handleTopologyParticleReaction(SCPUStateModel::topology_ref &topology,
                                                                   const SCPUEvaluateTopologyReactions::TREvent &event) {
    const auto& context = kernel->getKernelContext();
    const auto& top_registry = context.topology_registry();
    const auto& reaction = top_registry.spatial_reactions_by_type(event.t1, topology->type(), event.t2, topology_type_empty).at(event.reaction_idx);

    auto& model = kernel->getSCPUKernelStateModel();
    auto& data = *model.getParticleData();

    auto& entry1 = data.entry_at(event.idx1);
    auto& entry2 = data.entry_at(event.idx2);
    auto& entry1Type = entry1.type;
    auto& entry2Type = entry2.type;
    if(entry1Type == reaction.type1()) {
        entry1Type = reaction.type_to1();
        entry2Type = reaction.type_to2();
    } else {
        entry1Type = reaction.type_to2();
        entry2Type = reaction.type_to1();
    }
    if(event.topology_idx2 < 0) {
        entry1.topology_index = event.topology_idx;
        entry2.topology_index = event.topology_idx;

        if(reaction.is_fusion()) {
            topology->appendParticle(event.idx2, entry2Type, event.idx1, entry1Type);
        } else {
            topology->vertexForParticle(event.idx1)->setParticleType(entry1Type);
        }
    } else {
        throw std::logic_error("this branch should never be reached as topology-topology reactions are "
                                       "handeled in a different method");
    }
    topology->updateReactionRates(context.topology_registry().structural_reactions_of(topology->type()));
    topology->configure();
}

void SCPUEvaluateTopologyReactions::handleStructuralReaction(SCPUStateModel::topologies_vec &topologies,
                                                             std::vector<SCPUStateModel::topology> &new_topologies,
                                                             const SCPUEvaluateTopologyReactions::TREvent &event,
                                                             SCPUStateModel::topology_ref &topology) const {
    const auto &context = kernel->getKernelContext();
    const auto &reactions = context.topology_registry().structural_reactions_of(topology->type());
    const auto &reaction = reactions.at(static_cast<std::size_t>(event.reaction_idx));
    auto result = reaction.execute(*topology, kernel);
    if (!result.empty()) {
        // we had a topology fission, so we need to actually remove the current topology from the
        // data structure
        topologies.erase(topologies.begin() + event.topology_idx);
        //log::error("erased topology with index {}", event.topology_idx);
        assert(topology->isDeactivated());
        move(result.begin(), result.end(), back_inserter(new_topologies));
    } else {
        if (topology->isNormalParticle(*kernel)) {
            kernel->getSCPUKernelStateModel().getParticleData()->entry_at(topology->getParticles().front()).topology_index = -1;
            topologies.erase(topologies.begin() + event.topology_idx);
            //log::error("erased topology with index {}", event.topology_idx);
            assert(topology->isDeactivated());
        }
    }
}


}
}
}
}
}