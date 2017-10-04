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
 * @date 20.10.16
 */

#include <readdy/kernel/cpu/actions/reactions/ReactionUtils.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace actions {
namespace reactions {

data_t::DataUpdate handleEventsGillespie(
        CPUKernel *const kernel, scalar timeStep, bool filterEventsInAdvance, bool approximateRate,
        std::vector<event_t> &&events, std::vector<record_t> *maybeRecords, reaction_counts_t *maybeCounts) {
    using rdy_particle_t = readdy::model::Particle;
    const auto& fixPos = kernel->context().fixPositionFun();

    data_t::EntriesUpdate newParticles{};
    std::vector<data_t::size_type> decayedEntries {};

    if(!events.empty()) {
        const auto &ctx = kernel->context();
        auto data = kernel->getCPUKernelStateModel().getParticleData();
        /**
         * Handle gathered reaction events
         */
        {
            std::size_t nDeactivated = 0;
            const std::size_t nEvents = events.size();
            while (nDeactivated < nEvents) {
                const auto alpha = (*(events.end() - nDeactivated - 1)).cumulativeRate;
                const auto x = readdy::model::rnd::uniform_real<scalar>(static_cast<scalar>(0.), alpha);
                const auto eventIt = std::lower_bound(
                        events.begin(), events.end() - nDeactivated, x,
                        [](const event_t &elem1, scalar elem2) {
                            return elem1.cumulativeRate < elem2;
                        }
                );
                const auto event = *eventIt;
                if (eventIt == events.end() - nDeactivated) {
                    throw std::runtime_error("this should not happen (event not found)");
                }
                if (filterEventsInAdvance || shouldPerformEvent(event.reactionRate, timeStep, approximateRate)) {
                    /**
                     * Perform reaction
                     */
                    {

                        auto entry1 = event.idx1;
                        if (event.nEducts == 1) {
                            auto reaction = ctx.reactions().order1ByType(event.t1)[event.reactionIdx];
                            if(maybeRecords != nullptr) {
                                record_t record;
                                record.reactionIndex = event.reactionIdx;
                                performReaction(data, ctx, entry1, entry1, newParticles, decayedEntries, reaction, &record);
                                fixPos(record.where);
                                maybeRecords->push_back(record);
                            } else {
                                performReaction(data, ctx, entry1, entry1, newParticles, decayedEntries, reaction, nullptr);
                            }
                            if(maybeCounts != nullptr) {
                                auto &countsOrder1 = std::get<0>(*maybeCounts);
                                countsOrder1.at(event.t1).at(event.reactionIdx)++;
                            }
                        } else {
                            auto reaction = ctx.reactions().order2ByType(event.t1, event.t2)[event.reactionIdx];
                            if(maybeRecords != nullptr) {
                                record_t record;
                                record.reactionIndex = event.reactionIdx;
                                performReaction(data, ctx, entry1, event.idx2, newParticles, decayedEntries, reaction,
                                                &record);
                                fixPos(record.where);
                                maybeRecords->push_back(record);
                            } else {
                                performReaction(data, ctx, entry1, event.idx2, newParticles, decayedEntries, reaction,
                                                nullptr);
                            }
                            if(maybeCounts != nullptr) {
                                auto &countsOrder2 = std::get<1>(*maybeCounts);
                                countsOrder2.at(std::tie(event.t1, event.t2)).at(event.reactionIdx)++;
                            }
                        }
                    }
                    /**
                     * deactivate events whose educts have disappeared (including the just handled one)
                     */
                    {
                        auto _it = events.begin();
                        scalar cumsum = 0.0;
                        const auto idx1 = event.idx1;
                        if (event.nEducts == 1) {
                            while (_it < events.end() - nDeactivated) {
                                if ((*_it).idx1 == idx1 ||
                                    ((*_it).nEducts == 2 && (*_it).idx2 == idx1)) {
                                    ++nDeactivated;
                                    std::iter_swap(_it, events.end() - nDeactivated);
                                } else {
                                    cumsum += (*_it).reactionRate;
                                    (*_it).cumulativeRate = cumsum;
                                    ++_it;
                                }
                            }
                        } else {
                            const auto idx2 = event.idx2;
                            while (_it < events.end() - nDeactivated) {
                                if ((*_it).idx1 == idx1 || (*_it).idx1 == idx2 ||
                                    ((*_it).nEducts == 2 &&
                                     ((*_it).idx2 == idx1 || (*_it).idx2 == idx2))) {
                                    ++nDeactivated;
                                    std::iter_swap(_it, events.end() - nDeactivated);
                                } else {
                                    (*_it).cumulativeRate = cumsum;
                                    cumsum += (*_it).reactionRate;
                                    ++_it;
                                }
                            }
                        }

                    }
                } else {
                    nDeactivated++;
                    std::iter_swap(eventIt, events.end() - nDeactivated);
                    (*eventIt).cumulativeRate = (*eventIt).reactionRate;
                    if (eventIt > events.begin()) {
                        (*eventIt).cumulativeRate += (*(eventIt - 1)).cumulativeRate;
                    }
                    auto cumsum = (*eventIt).cumulativeRate;
                    for (auto _it = eventIt + 1; _it < events.end() - nDeactivated; ++_it) {
                        cumsum += (*_it).reactionRate;
                        (*_it).cumulativeRate = cumsum;
                    }
                }
            }
        }
    }
    return std::make_pair(std::move(newParticles), std::move(decayedEntries));
}
}
}
}
}
}