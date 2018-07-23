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
 * @file ReactionUtils.h<
 * @brief << brief description >>
 * @author clonker
 * @date 20.10.16
 */

#pragma once
#include <cmath>
#include <readdy/model/RandomProvider.h>
#include <readdy/kernel/cpu/CPUKernel.h>
#include <readdy/common/logging.h>
#include <readdy/common/boundary_condition_operations.h>
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
using reaction_counts_map = CPUStateModel::reaction_counts_map;
using neighbor_list = CPUStateModel::neighbor_list;

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
        std::vector<event_t> &&events, std::vector<record_t> *maybeRecords, reaction_counts_map *maybeCounts);

template<typename ParticleIndexCollection>
void gatherEvents(CPUKernel *const kernel, const ParticleIndexCollection &particles, const neighbor_list* nl,
                  const data_t *data, readdy::scalar &alpha, std::vector<event_t> &events) {
    const auto &box = kernel->context().boxSize();
    const auto &pbc = kernel->context().periodicBoundaryConditions();
    const auto& reaction_registry = kernel->context().reactions();
    for (const auto index : particles) {
        const auto &entry = data->entry_at(index);
        // this being false should really not happen, though
        if (!entry.deactivated) {
            // order 1
            {
                const auto &reactions = kernel->context().reactions().order1ByType(entry.type);
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
    for(std::size_t cell = 0; cell < nl->nCells(); ++cell) {
        for(auto particleIt = nl->particlesBegin(cell); particleIt != nl->particlesEnd(cell); ++particleIt) {
            const auto &idx1 = *particleIt;
            const auto &entry = data->entry_at(idx1);
            if(entry.deactivated) {
                log::critical("deactivated entry in neighbor list!");
                continue;
            }
            nl->forEachNeighbor(*particleIt, cell, [&](auto idx2) {
                if(idx1 > idx2) return;
                const auto &neighbor = data->entry_at(idx2);
                if(!neighbor.deactivated) {
                    const auto &reactions = kernel->context().reactions().order2ByType(entry.type, neighbor.type);
                    if (!reactions.empty()) {
                        const auto distSquared = bcs::distSquared(neighbor.pos, entry.pos, box, pbc);
                        for (auto itReactions = reactions.begin(); itReactions < reactions.end(); ++itReactions) {
                            const auto &react = *itReactions;
                            const auto rate = react->rate();
                            if (rate > 0 && distSquared < react->eductDistanceSquared()) {
                                alpha += rate;
                                events.emplace_back(2, react->nProducts(), idx1, idx2, rate, alpha,
                                                    static_cast<event_t::reaction_index_type>(itReactions -
                                                                                              reactions.begin()),
                                                    entry.type, neighbor.type);
                            }
                        }
                    }
                } else {
                    log::critical("deactivated entry in neighbor list!");
                }
            });
        }
    }
}

template<typename Reaction>
void performReaction(data_t* data, const readdy::model::Context& context, data_t::size_type idx1, data_t::size_type idx2,
                     data_t::EntriesUpdate& newEntries, std::vector<data_t::size_type>& decayedEntries,
                     Reaction* reaction, record_t* record) {

    const auto &pbc = context.periodicBoundaryConditions().data();
    const auto &box = context.boxSize().data();

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
            newEntries.emplace_back(bcs::applyPBC(entry1.pos - reaction->weight2() * reaction->productDistance() * n3,
                                                  box, pbc),
                                    reaction->products()[1], id);

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
            const auto difference = bcs::shortestDifference(e1Pos, e2Pos, box, pbc);
            if (reaction->educts()[0] == entry1.type) {
                newEntries.emplace_back(bcs::applyPBC(entry1.pos + reaction->weight1() * difference, box, pbc),
                                        reaction->products()[0], readdy::model::Particle::nextId());
            } else {
                newEntries.emplace_back(bcs::applyPBC(entry1.pos + reaction->weight2() * difference, box, pbc),
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
