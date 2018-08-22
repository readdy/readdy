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
 * @file CPUEvaluateTopologyReactions.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 15.06.17
 * @copyright GPL-3
 */

#include <readdy/kernel/cpu/actions/CPUEvaluateTopologyReactions.h>
#include <readdy/common/algorithm.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace actions {
namespace top {

CPUEvaluateTopologyReactions::CPUEvaluateTopologyReactions(CPUKernel *const kernel, scalar timeStep)
        : EvaluateTopologyReactions(timeStep), kernel(kernel) {}

/**
 * Struct holding information about a topology reaction event.
 */
struct CPUEvaluateTopologyReactions::TREvent {
    using index_type = CPUStateModel::data_type::size_type;

    rate_t cumulativeRate{0};
    rate_t rate{0};
    std::size_t topology_idx{0};
    // for topology-topology fusion only
    std::ptrdiff_t topology_idx2{-1};

    std::size_t reaction_idx{0};
    ParticleTypeId t1{0}, t2{0};
    // idx1 is always the particle that belongs to a topology
    index_type idx1{0}, idx2{0};
    bool spatial {false};

};

template<bool approximated>
bool performReactionEvent(scalar rate, scalar timeStep);

template<>
bool performReactionEvent<true>(const scalar rate, const scalar timeStep) {
    return readdy::model::rnd::uniform_real() < rate * timeStep;
}

template<>
bool performReactionEvent<false>(const scalar rate, const scalar timeStep) {
    return readdy::model::rnd::uniform_real() < 1 - std::exp(-rate * timeStep);
}

void CPUEvaluateTopologyReactions::perform(const util::PerformanceNode &node) {
    auto t = node.timeit();
    auto &model = kernel->getCPUKernelStateModel();
    const auto &context = kernel->context();
    auto &topologies = model.topologies();

    if (!topologies.empty()) {

        auto events = gatherEvents();

        if (!events.empty()) {

            std::vector<readdy::model::top::GraphTopology> new_topologies;

            {
                auto shouldEval = [this](const TREvent &event) {
                    return performReactionEvent<false>(event.rate, timeStep());
                };
                auto depending = [this](const TREvent &e1, const TREvent &e2) {
                    return eventsDependent(e1, e2);
                };
                auto eval = [&](const TREvent &event) {
                    auto &topology = topologies.at(event.topology_idx);
                    if (topology->isDeactivated()) {
                        throw std::logic_error(
                                fmt::format("deactivated topology with idx {} for {} event", event.topology_idx,
                                            event.spatial ? "spatial" : "structural"));
                    }
                    assert(!topology->isDeactivated());
                    if (!event.spatial) {
                        handleStructuralReaction(topologies, new_topologies, event, topology);
                    } else {
                        if (event.topology_idx2 >= 0) {
                            auto &top2 = topologies.at(static_cast<std::size_t>(event.topology_idx2));
                            if (topologyDeactivated(event.topology_idx2)) {
                                throw std::logic_error(
                                        fmt::format("encountered deactivated topology here, oh no (ix {}).",
                                                    event.topology_idx2));
                            }
                            handleTopologyTopologyReaction(topology, top2, event);
                        } else {
                            handleTopologyParticleReaction(topology, event);
                        }
                    }
                };

                /*auto postPerform = [&](const TREvent &event, std::size_t nDeactivated) {
                    for(auto it = events.begin(); it < events.end() - nDeactivated; ++it) {
                        if (eventsDependent(event, *it)) {
                            //log::critical("here ya go ye filth: {}", sss.str());
                            std::stringstream ss;
                            for(auto _it = events.begin(); _it < events.end() - nDeactivated; ++_it) {
                                ss << fmt::format("\tEvent: {}({}) + {}({})", _it->topology_idx, _it->idx1, _it->topology_idx2, _it->idx2) << "\n";
                            }
                            throw std::logic_error(fmt::format("events list contained event that shouldve been deactivated (t11={}, t12={}, t21={}, t22={}), (ix11={}, ix12={}, ix21={}, ix22={})\nlist:\n{}",
                                                               event.topology_idx, event.topology_idx2, it->topology_idx, it->topology_idx2, event.idx1, event.idx2, it->idx1, it->idx2, ss.str()));
                        }
                    }
                };*/
                algo::performEvents(events, shouldEval, depending, eval);
            }

            if (!new_topologies.empty()) {
                for (auto &&top : new_topologies) {
                    if (!top.isNormalParticle(*kernel)) {
                        // we have a new topology here, update data accordingly.
                        top.updateReactionRates(context.topologyRegistry().structuralReactionsOf(top.type()));
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

void CPUEvaluateTopologyReactions::handleStructuralReaction(CPUStateModel::topologies_vec &topologies,
                                                            std::vector<CPUStateModel::topology> &new_topologies,
                                                            const CPUEvaluateTopologyReactions::TREvent &event,
                                                            CPUStateModel::topology_ref &topology) const {
    const auto &topology_type_registry = kernel->context().topologyRegistry();
    auto &reaction = topology_type_registry.structuralReactionsOf(topology->type()).at(static_cast<std::size_t>(event.reaction_idx));
    auto result = reaction.execute(*topology, kernel);
    if (!result.empty()) {
        // we had a topology fission, so we need to actually remove the current topology from the
        // data structure
        topologies.erase(topologies.begin() + event.topology_idx);
        assert(topology->isDeactivated());

        for (auto &it : result) {
            if(!it.isNormalParticle(*kernel)) {
                new_topologies.push_back(std::move(it));
            }
        }
    } else {
        if (topology->isNormalParticle(*kernel)) {
            kernel->getCPUKernelStateModel().getParticleData()->entry_at(topology->getParticles().front()).topology_index = -1;
            topologies.erase(topologies.begin() + event.topology_idx);
            //log::error("erased topology with index {}", event.topology_idx);
            assert(topology->isDeactivated());
        }
    }
}

CPUEvaluateTopologyReactions::topology_reaction_events CPUEvaluateTopologyReactions::gatherEvents() {
    topology_reaction_events events;
    const auto &topology_types = kernel->context().topologyRegistry();
    {
        rate_t current_cumulative_rate = 0;
        std::size_t topology_idx = 0;
        for (auto &top : kernel->getCPUKernelStateModel().topologies()) {
            if (!top->isDeactivated()) {
                std::size_t reaction_idx = 0;
                for (const auto &reaction : topology_types.structuralReactionsOf(top->type())) {
                    TREvent event{};
                    event.rate = top->rates().at(reaction_idx);
                    event.cumulativeRate = event.rate + current_cumulative_rate;
                    current_cumulative_rate = event.cumulativeRate;
                    event.topology_idx = topology_idx;
                    event.reaction_idx = reaction_idx;

                    events.push_back(event);
                    ++reaction_idx;
                }
            }
            ++topology_idx;
        }

        const auto &context = kernel->context();

        static const CPUStateModel::topology_ref EMPTY_TOP {};

        if (!context.topologyRegistry().spatialReactionRegistry().empty()) {
            const auto &model = kernel->getCPUKernelStateModel();
            const auto &top_registry = context.topologyRegistry();
            const auto &box = context.boxSize().data();
            const auto &pbc = context.periodicBoundaryConditions().data();
            const auto &data = *kernel->getCPUKernelStateModel().getParticleData();
            const auto &nl = *kernel->getCPUKernelStateModel().getNeighborList();
            const auto &topologies = kernel->getCPUKernelStateModel().topologies();

            for (std::size_t cell = 0; cell < nl.nCells(); ++cell) {
                for(auto itParticle = nl.particlesBegin(cell); itParticle != nl.particlesEnd(cell); ++itParticle) {
                    const auto &entry = data.entry_at(*itParticle);
                    if (!entry.deactivated && top_registry.isSpatialReactionType(entry.type)) {
                        const auto entryTopologyDeactivated = topologyDeactivated(entry.topology_index);
                        const auto hasEntryTop = entry.topology_index >= 0 && !entryTopologyDeactivated;

                        nl.forEachNeighbor(*itParticle, cell, [&](std::size_t neighborIndex) {
                            const auto &neighbor = data.entry_at(neighborIndex);
                            const auto neighborTopDeactivated = topologyDeactivated(neighbor.topology_index);
                            const auto hasNeighborTop = neighbor.topology_index >= 0 && !neighborTopDeactivated;
                            if ((!hasEntryTop && !hasNeighborTop) || (hasNeighborTop && *itParticle > neighborIndex)) {
                                // use symmetry or skip entirely
                                return;
                            }
                            TopologyTypeId tt1 = hasEntryTop ? topologies.at(
                                    static_cast<std::size_t>(entry.topology_index))->type()
                                                                 : static_cast<TopologyTypeId>(-1);
                            TopologyTypeId tt2 = hasNeighborTop ? topologies.at(
                                    static_cast<std::size_t>(neighbor.topology_index))->type()
                                                                    : static_cast<TopologyTypeId>(-1);

                            const auto distSquared = bcs::distSquared(entry.pos, neighbor.pos, box, pbc);
                            std::size_t reaction_index = 0;
                            const auto &otherTop = hasNeighborTop ? model.topologies().at(
                                    static_cast<std::size_t>(neighbor.topology_index)
                            ) : EMPTY_TOP;
                            const auto &reactions = top_registry.spatialReactionsByType(entry.type, tt1,
                                                                                        neighbor.type, tt2);
                            for (const auto &reaction : reactions) {
                                if (!reaction.allow_self_connection() &&
                                    entry.topology_index == neighbor.topology_index) {
                                    ++reaction_index;
                                    continue;
                                }
                                if (distSquared < reaction.radius() * reaction.radius()) {
                                    TREvent event{};
                                    event.rate = reaction.rate();
                                    event.cumulativeRate = event.rate + current_cumulative_rate;
                                    current_cumulative_rate = event.cumulativeRate;
                                    if (hasEntryTop && !hasNeighborTop) {
                                        // entry is a topology, neighbor an ordinary particle
                                        event.topology_idx = static_cast<std::size_t>(entry.topology_index);
                                        event.t1 = entry.type;
                                        event.t2 = neighbor.type;
                                        event.idx1 = *itParticle;
                                        event.idx2 = neighborIndex;
                                    } else if (!hasEntryTop && hasNeighborTop) {
                                        // neighbor is a topology, entry an ordinary particle
                                        event.topology_idx = static_cast<std::size_t>(neighbor.topology_index);
                                        event.t1 = neighbor.type;
                                        event.t2 = entry.type;
                                        event.idx1 = neighborIndex;
                                        event.idx2 = *itParticle;
                                    } else if (hasEntryTop && hasNeighborTop) {
                                        // this is a topology-topology fusion
                                        event.topology_idx = static_cast<std::size_t>(entry.topology_index);
                                        event.topology_idx2 = static_cast<std::size_t>(neighbor.topology_index);
                                        event.t1 = entry.type;
                                        event.t2 = neighbor.type;
                                        event.idx1 = *itParticle;
                                        event.idx2 = neighborIndex;
                                    } else {
                                        log::critical("got no topology for topology-fusion");
                                    }
                                    event.reaction_idx = reaction_index;
                                    event.spatial = true;

                                    events.push_back(event);
                                }
                                ++reaction_index;
                            }
                        });
                    }
                }
            }
        }
    }
    return events;
}

void CPUEvaluateTopologyReactions::handleTopologyParticleReaction(CPUStateModel::topology_ref &topology,
                                                                  const CPUEvaluateTopologyReactions::TREvent &event) {
    const auto& context = kernel->context();
    const auto& top_registry = context.topologyRegistry();
    const auto& reaction = top_registry.spatialReactionsByType(event.t1, topology->type(), event.t2,
                                                               EmptyTopologyId).at(event.reaction_idx);

    auto& model = kernel->getCPUKernelStateModel();
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
            topology->vertexForParticle(event.idx1)->particleType() = entry1Type;
        }
    } else {
        throw std::logic_error("this branch should never be reached as topology-topology reactions are "
                                       "handeled in a different method");
    }
    if(topology->type() == reaction.top_type1()) {
        topology->type() = reaction.top_type_to1();
    } else {
        topology->type() = reaction.top_type_to2();
    }
    topology->updateReactionRates(context.topologyRegistry().structuralReactionsOf(topology->type()));
    topology->configure();
}

void CPUEvaluateTopologyReactions::handleTopologyTopologyReaction(CPUStateModel::topology_ref &t1,
                                                                  CPUStateModel::topology_ref &t2,
                                                                  const TREvent &event) {
    const auto& context = kernel->context();
    const auto& top_registry = context.topologyRegistry();
    const auto& reaction = top_registry.spatialReactionsByType(event.t1, t1->type(), event.t2, t2->type()).at(event.reaction_idx);

    auto& model = kernel->getCPUKernelStateModel();
    auto& data = *model.getParticleData();

    auto &entry1 = data.entry_at(event.idx1);
    auto &entry2 = data.entry_at(event.idx2);
    auto &entry1Type = entry1.type;
    auto &entry2Type = entry2.type;

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
            if(entry1Type == reaction.type1() && t1->type() == reaction.top_type1()) {
                v1->particleType() = reaction.type_to1();
                v2->particleType() = reaction.type_to2();
            } else {
                v1->particleType() = reaction.type_to2();
                v2->particleType() = reaction.type_to1();
            }
            if(!t1->graph().containsEdge(v1, v2)) {
                t1->graph().addEdge(v1, v2);
            }
            t1->type() = reaction.top_type_to1();
        } else {
            // merge topologies
            for(auto pidx : t2->getParticles()) {
                data.entry_at(pidx).topology_index = event.topology_idx;
            }
            auto &topologies = kernel->getCPUKernelStateModel().topologies();
            t1->appendTopology(*t2, event.idx2, entry2Type, event.idx1, entry1Type, reaction.top_type_to1());
            topologies.erase(topologies.begin() + event.topology_idx2);
        }
    } else {
        t1->vertexForParticle(event.idx1)->particleType() = entry1Type;
        t2->vertexForParticle(event.idx2)->particleType() = entry2Type;
        t1->type() = top_type_to1;
        t2->type() = top_type_to2;

        t2->updateReactionRates(context.topologyRegistry().structuralReactionsOf(t2->type()));
        t2->configure();
    }
    t1->updateReactionRates(context.topologyRegistry().structuralReactionsOf(t1->type()));
    t1->configure();
}

bool CPUEvaluateTopologyReactions::eventsDependent(const CPUEvaluateTopologyReactions::TREvent &evt1,
                                                   const CPUEvaluateTopologyReactions::TREvent &evt2) const {
    if(evt1.topology_idx == evt2.topology_idx || evt1.topology_idx == evt2.topology_idx2) {
        return true;
    }
    return evt1.topology_idx2 >= 0 && (evt1.topology_idx2 == evt2.topology_idx || evt1.topology_idx2 == evt2.topology_idx2);
}

bool CPUEvaluateTopologyReactions::topologyDeactivated(std::ptrdiff_t index) const {
    return index < 0 || kernel->getCPUKernelStateModel().topologies().at(
            static_cast<std::size_t>(index))->isDeactivated();
}

}
}
}
}
}
