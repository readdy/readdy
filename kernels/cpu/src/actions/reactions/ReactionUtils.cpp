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
 * @file ReactionUtils.cpp
 * @brief << brief description >>
 * @author clonker
 * @author chrisfroe
 * @date 20.10.16
 */

#include <readdy/kernel/cpu/actions/reactions/ReactionUtils.h>
#include <readdy/common/algorithm.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace actions {
namespace reactions {

data_t::DataUpdate handleEventsGillespie(
        CPUKernel *const kernel, scalar timeStep, bool filterEventsInAdvance, bool approximateRate,
        std::vector<event_t> &&events, std::vector<record_t> *maybeRecords, reaction_counts_map *maybeCounts) {
    using rdy_particle_t = readdy::model::Particle;
    const auto &box = kernel->context().boxSize().data();
    const auto &pbc = kernel->context().periodicBoundaryConditions().data();

    data_t::EntriesUpdate newParticles{};
    std::vector<data_t::size_type> decayedEntries{};

    if (!events.empty()) {
        const auto &ctx = kernel->context();
        auto data = kernel->getCPUKernelStateModel().getParticleData();
        /**
         * Handle gathered reaction events
         */
        {

            auto shouldEval = [&](const event_t &event) {
                return filterEventsInAdvance || shouldPerformEvent(event.rate, timeStep, approximateRate);
            };

            auto depending = [&](const event_t &e1, const event_t &e2) {
                return (e1.idx1 == e2.idx1 || (e2.nEducts == 2 && e1.idx1 == e2.idx2)
                        || (e1.nEducts == 2 && (e1.idx2 == e2.idx1 || (e2.nEducts == 2 && e1.idx2 == e2.idx2))));
            };

            auto eval = [&](const event_t &event) {
                auto entry1 = event.idx1;
                if (event.nEducts == 1) {
                    auto reaction = ctx.reactions().order1ByType(event.t1)[event.reactionIndex];
                    if (maybeRecords != nullptr) {
                        record_t record;
                        record.id = reaction->id();
                        performReaction(data, ctx, entry1, entry1, newParticles, decayedEntries, reaction, &record);
                        bcs::fixPosition(record.where, box, pbc);
                        maybeRecords->push_back(record);
                    } else {
                        performReaction(data, ctx, entry1, entry1, newParticles, decayedEntries, reaction, nullptr);
                    }
                    if (maybeCounts != nullptr) {
                        auto &counts = *maybeCounts;
                        counts.at(reaction->id())++;
                    }
                } else {
                    auto reaction = ctx.reactions().order2ByType(event.t1, event.t2)[event.reactionIndex];
                    if (maybeRecords != nullptr) {
                        record_t record;
                        record.id = reaction->id();
                        performReaction(data, ctx, entry1, event.idx2, newParticles, decayedEntries, reaction,
                                        &record);
                        bcs::fixPosition(record.where, box, pbc);
                        maybeRecords->push_back(record);
                    } else {
                        performReaction(data, ctx, entry1, event.idx2, newParticles, decayedEntries, reaction,
                                        nullptr);
                    }
                    if (maybeCounts != nullptr) {
                        auto &counts = *maybeCounts;
                        counts.at(reaction->id())++;
                    }
                }
            };

            algo::performEvents(events, shouldEval, depending, eval);
        }
    }
    return std::make_pair(std::move(newParticles), std::move(decayedEntries));
}
}
}
}
}
}
