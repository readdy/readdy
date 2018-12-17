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
 * @file algorithm.h
 * @brief << brief description >>
 * @author clonker
 * @date 3/1/18
 */


#pragma once

#include <vector>
#include <algorithm>

#include "common.h"
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
                                 const EvaluateOnTopology &evaluateOnTopology) {
    // Evaluate on particles
    {
        std::for_each(particleContainer.begin(), particleContainer.end(), [&](auto &&entry){
            if (!entry.deactivated) {
                evaluateOnParticle(entry);
            }
        });
    }

    // Evaluate on interactions
    {
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
        for (auto &&topology : topologyContainer) {
            if (!topology->isDeactivated()) {
                evaluateOnTopology(topology);
            }
        }
    }
}


}
}
