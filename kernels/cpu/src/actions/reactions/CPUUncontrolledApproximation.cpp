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
 * @file UncontrolledApproximation.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 20.10.16
 */

#include <future>
#include <random>

#include <readdy/kernel/cpu/actions/reactions/CPUUncontrolledApproximation.h>
#include <readdy/kernel/cpu/actions/reactions/Event.h>
#include <readdy/kernel/cpu/actions/reactions/ReactionUtils.h>
#include <readdy/kernel/cpu/data/DataContainer.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace actions {
namespace reactions {

using data_t = data::EntryDataContainer;
using event_t = Event;
using data_iter_t = data_t::const_iterator;
using neighbor_list = CPUStateModel::neighbor_list;
using nl_bounds = std::tuple<std::size_t, std::size_t>;
using entry_type = data_t::Entries::value_type;

using event_future_t = std::future<std::vector<event_t>>;
using event_promise_t = std::promise<std::vector<event_t>>;

CPUUncontrolledApproximation::CPUUncontrolledApproximation(CPUKernel *const kernel, scalar timeStep)
        : super(timeStep), kernel(kernel) {

}

void findEvents(std::size_t /*tid*/, data_iter_t begin, data_iter_t end, nl_bounds nlBounds,
                const CPUKernel *const kernel, scalar dt, bool approximateRate, const neighbor_list &nl,
                event_promise_t &events, std::promise<std::size_t> &n_events) {
    std::vector<event_t> eventsUpdate;
    const auto &data = *kernel->getCPUKernelStateModel().getParticleData();
    const auto &box = kernel->context().boxSize().data();
    const auto &pbc = kernel->context().periodicBoundaryConditions().data();
    auto index = static_cast<std::size_t>(std::distance(data.begin(), begin));
    for (auto it = begin; it != end; ++it, ++index) {
        const auto &entry = *it;
        // this being false should really not happen, though
        if (!entry.deactivated) {
            // order 1
            {
                const auto &reactions = kernel->context().reactions().order1ByType(entry.type);
                for (auto it_reactions = reactions.begin(); it_reactions != reactions.end(); ++it_reactions) {
                    const auto rate = (*it_reactions)->rate();
                    if (rate > 0 && shouldPerformEvent(rate, dt, approximateRate)) {
                        eventsUpdate.emplace_back(1, (*it_reactions)->nProducts(), index, index, rate, 0,
                                                  static_cast<event_t::reaction_index_type>(it_reactions -
                                                                                            reactions.begin()),
                                                  entry.type, 0);
                    }
                }
            }
        }
    }
    for(auto cell = std::get<0>(nlBounds); cell != std::get<1>(nlBounds); ++cell) {
        for(auto particleIt = nl.particlesBegin(cell); particleIt != nl.particlesEnd(cell); ++particleIt) {
            const auto &entry = data.entry_at(*particleIt);
            if(entry.deactivated) {
                log::critical("deactivated entry in uncontrolled approximation!");
                continue;
            }

            nl.forEachNeighbor(*particleIt, cell, [&](const auto neighborIdx) {
                const auto &neighbor = data.entry_at(neighborIdx);
                if(!neighbor.deactivated) {
                    const auto &reactions = kernel->context().reactions().order2ByType(entry.type, neighbor.type);
                    if (!reactions.empty()) {
                        const auto distSquared = bcs::distSquared(neighbor.pos, entry.pos, box, pbc);
                        for (auto it_reactions = reactions.begin(); it_reactions < reactions.end(); ++it_reactions) {
                            const auto &react = *it_reactions;
                            const auto rate = react->rate();
                            if (rate > 0 && distSquared < react->eductDistanceSquared()
                                && shouldPerformEvent(rate, dt, approximateRate)) {
                                const auto reaction_index = static_cast<event_t::reaction_index_type>(it_reactions -
                                                                                                      reactions.begin());
                                eventsUpdate.emplace_back(2, react->nProducts(), *particleIt, neighborIdx,
                                                          rate, 0, reaction_index, entry.type, neighbor.type);
                            }
                        }
                    }
                }
            });
        }
    }

    n_events.set_value(eventsUpdate.size());
    events.set_value(std::move(eventsUpdate));
}

void CPUUncontrolledApproximation::perform(const util::PerformanceNode &node) {
    auto t = node.timeit();
    const auto &ctx = kernel->context();
    auto &stateModel = kernel->getCPUKernelStateModel();
    auto nl = stateModel.getNeighborList();
    auto &data = nl->data();
    const auto &box = ctx.boxSize().data();
    const auto &pbc = ctx.periodicBoundaryConditions().data();

    if (ctx.recordReactionsWithPositions()) {
        kernel->getCPUKernelStateModel().reactionRecords().clear();
    }
    if (ctx.recordReactionCounts()) {
        stateModel.resetReactionCounts();
    }

    // gather events
    std::vector<std::promise<std::size_t>> n_events_promises(kernel->getNThreads());
    std::vector<event_promise_t> promises(kernel->getNThreads());
    {

        auto &pool = kernel->pool();

        std::vector<std::function<void(std::size_t)>> executables;
        executables.reserve(kernel->getNThreads());

        std::size_t grainSize = data.size() / kernel->getNThreads();
        std::size_t nlGrainSize = nl->nCells() / kernel->getNThreads();

        auto it = data.cbegin();
        std::size_t it_nl = 0;
        for (auto i = 0U; i < kernel->getNThreads()-1 ; ++i) {
            auto itNext = std::min(it+grainSize, data.cend());

            auto nlNext = std::min(it_nl + nlGrainSize, nl->nCells());
            auto bounds_nl = std::make_tuple(it_nl, nlNext);

            pool.push(findEvents, it, itNext, bounds_nl, kernel, timeStep(), false, std::cref(*nl),
                      std::ref(promises.at(i)), std::ref(n_events_promises.at(i)));

            it = itNext;
            it_nl = nlNext;
        }
        pool.push(findEvents, it, data.cend(), std::make_tuple(it_nl, nl->nCells()), kernel, timeStep(), false,
                  std::cref(*nl), std::ref(promises.back()), std::ref(n_events_promises.back()));
    }

    // collect events
    std::vector<event_t> events;
    {
        std::size_t n_events = 0;
        for (auto &&f : n_events_promises) {
            n_events += f.get_future().get();
        }
        events.reserve(n_events);
        for (auto &&f : promises) {
            auto eventUpdate = std::move(f.get_future().get());
            auto mBegin = std::make_move_iterator(eventUpdate.begin());
            auto mEnd = std::make_move_iterator(eventUpdate.end());
            events.insert(events.end(), mBegin, mEnd);
        }
    }

    // shuffle reactions
    std::shuffle(events.begin(), events.end(), std::mt19937(std::random_device()()));

    // execute reactions
    {
        data_t::EntriesUpdate newParticles{};
        std::vector<data_t::size_type> decayedEntries{};

        // todo better conflict detection?
        for (auto it = events.begin(); it != events.end(); ++it) {
            auto &event = *it;
            if (event.cumulativeRate == 0) {
                auto entry1 = event.idx1;
                if (event.nEducts == 1) {
                    auto reaction = ctx.reactions().order1ByType(event.t1)[event.reactionIndex];
                    if (ctx.recordReactionsWithPositions()) {
                        record_t record;
                        record.id = reaction->id();
                        performReaction(&data, ctx, entry1, entry1, newParticles, decayedEntries, reaction, &record);
                        bcs::fixPosition(record.where, box, pbc);
                        kernel->getCPUKernelStateModel().reactionRecords().push_back(record);
                    } else {
                        performReaction(&data, ctx, entry1, entry1, newParticles, decayedEntries, reaction, nullptr);
                    }
                    if (ctx.recordReactionCounts()) {
                        auto &counts = stateModel.reactionCounts();
                        counts.at(reaction->id())++;
                    }
                    for (auto _it2 = it + 1; _it2 != events.end(); ++_it2) {
                        if (_it2->idx1 == entry1 || _it2->idx2 == entry1) {
                            _it2->cumulativeRate = 1;
                        }
                    }
                } else {
                    auto reaction = ctx.reactions().order2ByType(event.t1, event.t2)[event.reactionIndex];
                    if (ctx.recordReactionsWithPositions()) {
                        record_t record;
                        record.id = reaction->id();
                        performReaction(&data, ctx, entry1, event.idx2, newParticles, decayedEntries, reaction, &record);
                        bcs::fixPosition(record.where, box, pbc);
                        kernel->getCPUKernelStateModel().reactionRecords().push_back(record);
                    } else {
                        performReaction(&data, ctx, entry1, event.idx2, newParticles, decayedEntries, reaction, nullptr);
                    }
                    if (ctx.recordReactionCounts()) {
                        auto &counts = stateModel.reactionCounts();
                        counts.at(reaction->id())++;
                    }
                    for (auto _it2 = it + 1; _it2 != events.end(); ++_it2) {
                        if (_it2->idx1 == entry1 || _it2->idx2 == entry1 ||
                            _it2->idx1 == event.idx2 || _it2->idx2 == event.idx2) {
                            _it2->cumulativeRate = 1;
                        }
                    }
                }
            }
        }
        data.update(std::make_pair(std::move(newParticles), std::move(decayedEntries)));
    }
}
}
}
}
}
}
