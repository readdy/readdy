/********************************************************************
 * Copyright © 2018 Computational Molecular Biology Group,          *
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
 * @file algorithm.h
 * @brief << brief description >>
 * @author clonker
 * @date 3/1/18
 */


#pragma once

#include <vector>
#include <algorithm>

#include "common.h"
#include <readdy/common/Timer.h>
#include "../model/RandomProvider.h"

namespace readdy {
namespace algo {

namespace detail {
template<typename Events>
void noPostPerform(const typename Events::value_type &, std::size_t) {}
}

template<typename Events, typename ShouldEvaluate, typename Depending, typename Evaluate,
        typename PostPerform = std::function<void(typename Events::value_type, std::size_t)>>
inline void performEvents(Events &events, const ShouldEvaluate &shouldEvaluate, const Depending &depending,
                   const Evaluate &evaluate, const PostPerform &postPerform = detail::noPostPerform<Events>) {
    if(!events.empty()) {
        std::size_t nDeactivated = 0;
        const std::size_t nEvents = events.size();

        {
            // update cumulative rates
            scalar cumulativeRate = c_::zero;
            for(auto &event : events) {
                event.cumulativeRate = event.rate + cumulativeRate;
                cumulativeRate = event.cumulativeRate;
            }
        }

        while(nDeactivated < nEvents) {
            const auto cumulativeRate = (std::end(events) - nDeactivated - 1)->cumulativeRate;
            const auto x = readdy::model::rnd::uniform_real(c_::zero, cumulativeRate);
            const auto eventIt = std::lower_bound(
                    std::begin(events), std::end(events) - nDeactivated, x, [](const auto &elem1, const auto elem2) {
                        return elem1.cumulativeRate < elem2;
                    }
            );

            if (eventIt != events.end()) {

                if(shouldEvaluate(*eventIt)) {
                    evaluate(*eventIt);

                    auto evaluatedEvent = *eventIt;

                    // shift all events to the end that depend on this particular one
                    auto _it = std::begin(events);
                    scalar cumsum = 0.0;
                    while (_it < std::end(events) - nDeactivated) {
                        if (depending(*_it, evaluatedEvent)) {
                            ++nDeactivated;
                            std::iter_swap(_it, events.end() - nDeactivated);
                        } else {
                            cumsum += _it->rate;
                            _it->cumulativeRate = cumsum;
                            ++_it;
                        }
                    }

                    postPerform(evaluatedEvent, nDeactivated);
                } else {
                    // remove event from the list (ie shift it to the end)
                    ++nDeactivated;

                    // swap with last element that is not yet deactivated
                    std::iter_swap(eventIt, events.end() - nDeactivated);

                    // this element's cumulative rate gets initialized with its own rate
                    eventIt->cumulativeRate = eventIt->rate;
                    if (eventIt > std::begin(events)) {
                        // and then increased by the cumulative rate of its predecessor
                        eventIt->cumulativeRate += (eventIt - 1)->cumulativeRate;
                    }
                    // now update the cumulative rates of all following elements...
                    auto cumsum = (*eventIt).cumulativeRate;
                    for (auto _it = eventIt + 1; _it < std::end(events) - nDeactivated; ++_it) {
                        cumsum += (*_it).rate;
                        (*_it).cumulativeRate = cumsum;
                    }
                }

            } else {
                throw std::logic_error("internal error: performEvents - no event could be selected");
            }
        }
    }
}

template<typename ParticleContainer, typename EvaluateOnParticle, typename InteractionContainer,
        typename EvaluateOnInteraction, typename TopologyContainer, typename EvaluateOnTopology>
inline void evaluateOnContainers(ParticleContainer &&particleContainer,
                                 const EvaluateOnParticle &evaluateOnParticle,
                                 InteractionContainer &&interactionContainer,
                                 const EvaluateOnInteraction &evaluateOnInteraction,
                                 TopologyContainer &&topologyContainer,
                                 const EvaluateOnTopology &evaluateOnTopology,
                                 const util::PerformanceNode &node) {
    // Evaluate on particles
    {
        auto tEvaluateOnParticles = node.subnode("evaluate on particles").timeit();
        std::for_each(particleContainer.begin(), particleContainer.end(), [&](auto &&entry){
            if (!entry.deactivated) {
                evaluateOnParticle(entry);
            }
        });
    }

    // Evaluate on interactions
    {
        auto tEvaluateOnInteractions = node.subnode("evaluate on interactions").timeit();
        for (auto cell = 0_z; cell < interactionContainer.nCells(); ++cell) {
            for (auto it = interactionContainer.particlesBegin(cell); it != interactionContainer.particlesEnd(cell); ++it) {
                auto pidx = *it;
                auto &&entry = particleContainer.entry_at(pidx);
                interactionContainer.forEachNeighbor(it, cell, [&](const std::size_t neighbor) {
                    auto &&neighborEntry = particleContainer.entry_at(neighbor);
                    evaluateOnInteraction(entry, neighborEntry);

                });
            }
        }
    }

    // Evaluate on topologies
    {
        auto tTopologies = node.subnode("evaluate on topologies").timeit();
        for (auto &&topology : topologyContainer) {
            if (!topology->isDeactivated()) {
                evaluateOnTopology(topology);
            }
        }
    }
}


}
}
