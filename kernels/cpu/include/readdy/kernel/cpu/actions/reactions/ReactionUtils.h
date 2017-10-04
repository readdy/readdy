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
 * @file ReactionUtils.h<
 * @brief << brief description >>
 * @author clonker
 * @date 20.10.16
 */

#pragma once
#include <cmath>
#include <readdy/model/RandomProvider.h>
#include <readdy/kernel/cpu/CPUKernel.h>
#include <readdy/kernel/cpu/nl/AdaptiveNeighborList.h>
#include <readdy/common/logging.h>
#include "Event.h"

namespace readdy {
namespace kernel {
namespace cpu {
namespace actions {
namespace reactions {

using data_t = data::EntryDataContainer;
using cpu_kernel = readdy::kernel::cpu::CPUKernel;
using reaction_type = readdy::model::reactions::ReactionType;
using ctx_t = std::remove_const<decltype(std::declval<cpu_kernel>().context())>::type;
using event_t = Event;
using record_t = readdy::model::reactions::ReactionRecord;
using reaction_counts_order1_map = readdy::model::observables::ReactionCounts::reaction_counts_order1_map;
using reaction_counts_order2_map = readdy::model::observables::ReactionCounts::reaction_counts_order2_map;
using reaction_counts_t = std::pair<reaction_counts_order1_map, reaction_counts_order2_map>;

template<bool approximated>
bool performReactionEvent(const readdy::scalar rate, const readdy::scalar timeStep) {
    if (approximated) {
        return readdy::model::rnd::uniform_real<scalar>() < rate * timeStep;
    } else {
        return readdy::model::rnd::uniform_real<scalar>() < 1 - std::exp(-rate * timeStep);
    }
}


inline bool shouldPerformEvent(const readdy::scalar rate, const readdy::scalar timestep, bool approximated) {
    return approximated ? performReactionEvent<true>(rate, timestep) : performReactionEvent<false>(rate, timestep);
}

data_t::DataUpdate handleEventsGillespie(
        CPUKernel* kernel, readdy::scalar timeStep,
        bool filterEventsInAdvance, bool approximateRate,
        std::vector<event_t> &&events, std::vector<record_t> *maybeRecords, reaction_counts_t *maybeCounts);

template<typename ParticleIndexCollection>
void gatherEvents(CPUKernel *const kernel, const ParticleIndexCollection &particles, const neighbor_list* nl,
                  const data_t *data, readdy::scalar &alpha, std::vector<event_t> &events,
                  const readdy::model::Context::dist_squared_fun& d2) {
    const auto& reaction_registry = kernel->context().reactions();
    for (const auto index : particles) {
        const auto &entry = data->entry_at(index);
        // this being false should really not happen, though
        if (!entry.deactivated) {
            // order 1
            {
                const auto &reactions = kernel->context().reactions().order1_by_type(entry.type);
                for (auto it = reactions.begin(); it != reactions.end(); ++it) {
                    const auto rate = (*it)->rate();
                    if (rate > 0) {
                        alpha += rate;
                        events.emplace_back(
                                1, (*it)->nProducts(), index, 0, rate, alpha,
                                static_cast<event_t::reaction_index_type>(it - reactions.begin()),
                                entry.type, 0);
                    }
                }
            }
        }
    }
    // order 2
    for (auto it : *nl) {
        auto index = it.current_particle();
        const auto &entry = data->entry_at(index);
        if(!entry.deactivated) {
            for (auto idx_neighbor : it) {
                if (index > idx_neighbor) continue;
                const auto &neighbor = data->entry_at(idx_neighbor);
                if (!neighbor.deactivated) {
                    const auto &reactions = kernel->context().reactions().order2_by_type(entry.type,
                                                                                                  neighbor.type);
                    if (!reactions.empty()) {
                        const auto distSquared = d2(neighbor.pos, entry.pos);
                        for (auto itReactions = reactions.begin(); itReactions < reactions.end(); ++itReactions) {
                            const auto &react = *itReactions;
                            const auto rate = react->rate();
                            if (rate > 0 && distSquared < react->eductDistanceSquared()) {
                                alpha += rate;
                                events.emplace_back(2, react->nProducts(), index, idx_neighbor, rate, alpha,
                                                    static_cast<event_t::reaction_index_type>(itReactions -
                                                                                              reactions.begin()),
                                                    entry.type, neighbor.type);
                            }
                        }
                    }
                } else {
                    log::critical("deactivated entry in neighbor list!");
                }
            }
        }
    }
}

template<typename Reaction>
void performReaction(data_t* data, const readdy::model::Context& context, data_t::size_type idx1, data_t::size_type idx2,
                     data_t::EntriesUpdate& newEntries, std::vector<data_t::size_type>& decayedEntries,
                     Reaction* reaction, record_t* record) {
    const auto& pbc = context.applyPBCFun();
    auto& entry1 = data->entry_at(idx1);
    auto& entry2 = data->entry_at(idx2);
    if(record) {
        record->type = static_cast<int>(reaction->type());
        record->where = (entry1.pos + entry2.pos) / 2.;
        record->educts[0] = entry1.id;
        record->educts[1] = entry2.id;
        record->types_from[0] = entry1.type;
        record->types_from[1] = entry2.type;
    }
    switch(reaction->type()) {
        case reaction_type::Decay: {
            decayedEntries.push_back(idx1);
            break;
        }
        case reaction_type::Conversion: {
            entry1.type = reaction->products()[0];
            entry1.id = readdy::model::Particle::nextId();
            if(record) record->products[0] = entry1.id;
            break;
        }
        case reaction_type::Enzymatic: {
            if (entry1.type == reaction->educts()[1]) {
                // p1 is the catalyst
                entry2.type = reaction->products()[0];
                entry2.id = readdy::model::Particle::nextId();
            } else {
                // p2 is the catalyst
                entry1.type = reaction->products()[0];
                entry1.id = readdy::model::Particle::nextId();
            }
            if(record) {
                record->products[0] = entry1.id;
                record->products[1] = entry2.id;
            }
            break;
        }
        case reaction_type::Fission: {
            auto n3 = readdy::model::rnd::normal3<readdy::scalar>(0, 1);
            n3 /= std::sqrt(n3 * n3);

            //readdy::model::Particle p (, reaction->products()[1]);
            const auto id = readdy::model::Particle::nextId();
            newEntries.emplace_back(pbc(entry1.pos - reaction->weight2() * reaction->productDistance() * n3), reaction->products()[1], id);

            entry1.type = reaction->products()[0];
            entry1.id = readdy::model::Particle::nextId();
            data->displace(idx1, reaction->weight1() * reaction->productDistance() * n3);
            if(record) {
                record->products[0] = entry1.id;
                record->products[1] = id;
            }
            break;
        }
        case reaction_type::Fusion: {
            const auto& e1Pos = entry1.pos;
            const auto& e2Pos = entry2.pos;
            if (reaction->educts()[0] == entry1.type) {
                newEntries.emplace_back(pbc(entry1.pos + reaction->weight1() * (e2Pos - e1Pos)),
                                        reaction->products()[0], readdy::model::Particle::nextId());
            } else {
                newEntries.emplace_back(pbc(entry1.pos + reaction->weight2() * (e2Pos - e1Pos)),
                                        reaction->products()[0], readdy::model::Particle::nextId());
            }
            decayedEntries.push_back(idx1);
            decayedEntries.push_back(idx2);
            if(record) record->products[0] = newEntries.back().id;
            break;
        }
    }
}

}
}
}
}
}
