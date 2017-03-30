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
 * @date 23.11.16
 */

#include <future>

#include <readdy/common/thread/scoped_thread.h>

#include <readdy/kernel/cpu_dense/actions/reactions/CPUDUncontrolledApproximation.h>
#include <readdy/kernel/cpu_dense/actions/reactions/Event.h>
#include <readdy/kernel/cpu_dense/actions/reactions/ReactionUtils.h>

namespace readdy {
namespace kernel {
namespace cpu_dense {
namespace actions {
namespace reactions {

namespace thd = readdy::util::thread;

using event_t = Event;
using neighbor_list_t = model::CPUDNeighborList;
using data_t = neighbor_list_t::data_t;
using data_iter_t = data_t::const_iterator;
using neighbor_list_iter_t = neighbor_list_t::const_iterator;
using entry_type = data_t::Entry;

using event_future_t = std::future<std::vector<event_t>>;
using event_promise_t = std::promise<std::vector<event_t>>;

CPUDUncontrolledApproximation::CPUDUncontrolledApproximation(const CPUDKernel *const kernel, double timeStep)
        : super(timeStep), kernel(kernel) {

}

void findEvents(data_iter_t begin, data_iter_t end, neighbor_list_iter_t nl_begin, const CPUDKernel *const kernel,
                double timeStep, bool approximateRate, event_promise_t events, std::promise<std::size_t> n_events) {
    std::vector<event_t> eventsUpdate;
    const auto &data = *kernel->getCPUDKernelStateModel().getParticleData();
    auto it = begin;
    auto it_nl = nl_begin;
    auto index = static_cast<std::size_t>(std::distance(data.begin(), begin));
    for (; it != end; ++it, ++it_nl, ++index) {
        const auto &entry = *it;
        // this being false should really not happen, though
        if (!entry.deactivated) {
            // order 1
            {
                const auto &reactions = kernel->getKernelContext().reactions().order1_by_type(entry.type);
                for (auto it_reactions = reactions.begin(); it_reactions != reactions.end(); ++it_reactions) {
                    const auto rate = (*it_reactions)->getRate();
                    if (rate > 0 && shouldPerformEvent(rate, timeStep, approximateRate)) {
                        eventsUpdate.push_back(
                                {1, (*it_reactions)->getNProducts(), index, index, rate, 0,
                                 static_cast<event_t::reaction_index_type>(it_reactions - reactions.begin()),
                                 entry.type, 0});
                    }
                }
            }
            // order 2
            for(auto it_neighbors = it_nl->cbegin(); it_neighbors != it_nl->cend(); ++it_neighbors) {
                const auto neighbor_index = it_neighbors->idx;
                if (index > neighbor_index) continue;
                const auto &neighbor = data.entry_at(neighbor_index);
                const auto &reactions = kernel->getKernelContext().reactions().order2_by_type(entry.type,
                                                                                                     neighbor.type);
                if (!reactions.empty()) {
                    const auto distSquared = it_neighbors->d2;
                    for (auto it_reactions = reactions.begin(); it_reactions < reactions.end(); ++it_reactions) {
                        const auto &react = *it_reactions;
                        const auto rate = react->getRate();
                        if (rate > 0 && distSquared < react->getEductDistanceSquared()
                            && shouldPerformEvent(rate, timeStep, approximateRate)) {
                            const auto reaction_index = static_cast<event_t::reaction_index_type>(it_reactions -
                                                                                                  reactions.begin());
                            eventsUpdate.push_back({2, react->getNProducts(), index, neighbor_index, rate, 0,
                                                    reaction_index, entry.type, neighbor.type});
                        }
                    }
                }
            }
        }
    }

    n_events.set_value(eventsUpdate.size());
    events.set_value(std::move(eventsUpdate));
}

void CPUDUncontrolledApproximation::perform() {
    const auto &ctx = kernel->getKernelContext();
    const auto &fixPos = ctx.getFixPositionFun();
    const auto &dt = timeStep;
    auto &data = *kernel->getCPUDKernelStateModel().getParticleData();
    auto &nl = *kernel->getCPUDKernelStateModel().getNeighborList();

    // gather events
    std::vector<std::future<std::size_t>> n_eventsFutures;
    std::vector<event_future_t> eventFutures;
    {
        const std::size_t grainSize = data.size() / kernel->getNThreads();

        std::vector<thd::scoped_thread> threads;

        auto it = data.cbegin();
        auto it_nl = nl.cbegin();
        for (unsigned int i = 0; i < kernel->getNThreads()-1; ++i) {
            event_promise_t eventPromise;
            eventFutures.push_back(eventPromise.get_future());
            std::promise<std::size_t> n_events;
            n_eventsFutures.push_back(n_events.get_future());

            threads.push_back(thd::scoped_thread(
                    std::thread(findEvents, it, it + grainSize, it_nl, kernel, timeStep, true,
                                std::move(eventPromise), std::move(n_events))
            ));
            std::advance(it, grainSize);
            std::advance(it_nl, grainSize);
        }
        {
            event_promise_t eventPromise;
            eventFutures.push_back(eventPromise.get_future());
            std::promise<std::size_t> n_events;
            n_eventsFutures.push_back(n_events.get_future());

            threads.push_back(thd::scoped_thread(
                    std::thread(findEvents, it, data.cend(), it_nl, kernel, timeStep, true,
                                std::move(eventPromise), std::move(n_events))
            ));
        }
    }

    // collect events
    std::vector<event_t> events;
    {
        std::size_t n_events = 0;
        for (auto &&f : n_eventsFutures) {
            n_events += f.get();
        }
        events.reserve(n_events);
        for (auto &&f : eventFutures) {
            auto eventUpdate = std::move(f.get());
            auto mBegin = std::make_move_iterator(eventUpdate.begin());
            auto mEnd = std::make_move_iterator(eventUpdate.end());
            events.insert(events.end(), mBegin, mEnd);
        }
    }

    // shuffle reactions
    std::random_shuffle(events.begin(), events.end());

    // execute reactions
    {
        data_t::update_t newParticles{};
        std::vector<data_t::index_t> decayedEntries {};

        for(auto it = events.begin(); it != events.end(); ++it) {
            auto& event = *it;
            if(event.cumulativeRate == 0) {
                auto entry1 = event.idx1;
                if (event.nEducts == 1) {
                    auto reaction = ctx.reactions().order1_by_type(event.t1)[event.reactionIdx];
                    performReaction<true>(data, entry1, entry1, newParticles, reaction);
                } else {
                    auto reaction = ctx.reactions().order2_by_type(event.t1, event.t2)[event.reactionIdx];
                    performReaction<true>(data, entry1, event.idx2, newParticles, reaction);
                }
            }
        }

        data.deactivateMarked();
        data.update(std::move(newParticles));
    }
}
}
}
}
}
}