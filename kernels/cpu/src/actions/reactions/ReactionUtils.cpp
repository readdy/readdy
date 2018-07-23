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
