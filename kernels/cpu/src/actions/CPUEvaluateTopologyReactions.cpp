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
 * @file CPUEvaluateTopologyReactions.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 15.06.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <readdy/kernel/cpu/actions/CPUEvaluateTopologyReactions.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace actions {
namespace top {


CPUEvaluateTopologyReactions::CPUEvaluateTopologyReactions(CPUKernel *const kernel, double timeStep)
        : EvaluateTopologyReactions(timeStep), kernel(kernel){}

template<bool approximated>
bool performReactionEvent(const double rate, const double timeStep) {
    if (approximated) {
        return readdy::model::rnd::uniform_real() < rate * timeStep;
    } else {
        return readdy::model::rnd::uniform_real() < 1 - std::exp(-rate * timeStep);
    }
}

bool shouldPerformEvent(const double rate, const double timestep, bool approximated) {
    return approximated ? performReactionEvent<true>(rate, timestep) : performReactionEvent<false>(rate, timestep);
}

void CPUEvaluateTopologyReactions::perform() {
    using rate_t = readdy::model::top::GraphTopology::rate_t;
    auto& topologies = kernel->getCPUKernelStateModel().topologies();

    struct TREvent {
        rate_t cumulative_rate;
        rate_t own_rate;
        std::size_t topology_idx;
        std::size_t reaction_idx;
    };

    std::vector<TREvent> events;

    {
        rate_t current_cumulative_rate = 0;
        std::size_t topology_idx = 0;
        for (auto &topRef : topologies) {

            if (!topRef->isDeactivated()) {

                std::size_t reaction_idx = 0;
                for (const auto &reaction : topRef->registeredReactions()) {

                    TREvent event;
                    event.own_rate = std::get<1>(reaction);
                    event.cumulative_rate = event.own_rate + current_cumulative_rate;
                    current_cumulative_rate = event.cumulative_rate;
                    event.topology_idx = topology_idx;
                    event.reaction_idx = reaction_idx;

                    events.push_back(std::move(event));

                    ++reaction_idx;
                }

            }
            ++topology_idx;
        }
    }

    if(!events.empty()) {

        std::size_t n_processed = 0;
        const std::size_t n_events = events.size();

        std::vector<readdy::model::top::GraphTopology> new_topologies;
        while(n_processed < n_events) {
            const auto cumulative_rate = events.at(n_events - n_processed - 1).cumulative_rate;

            const auto x = readdy::model::rnd::uniform_real(0., cumulative_rate);

            const auto eventIt = std::lower_bound(
                    events.begin(), events.end() - n_processed, x, [](const TREvent &elem1, rate_t elem2) {
                        return elem1.cumulative_rate < elem2;
                    }
            );

            if(eventIt != events.end()) {
                const auto& event = *eventIt;

                if(shouldPerformEvent(event.own_rate, timeStep, true)) {
                    // perform the event!
                    auto& topology = topologies.at(event.topology_idx);
                    assert(!topology->isDeactivated());
                    auto& reaction = topology->registeredReactions().at(event.reaction_idx);
                    auto result = std::get<0>(reaction).execute(*topology, kernel);
                    if(!result.empty()) {
                        // we had a topology fission, so we need to actually remove the current topology from the
                        // data structure
                        topologies.erase(topologies.begin() + event.topology_idx);
                        assert(topology->isDeactivated());
                        std::move(result.begin(), result.end(), std::back_inserter(new_topologies));
                    }
                    // as in the non-performing case, deactivate this event and move it to the back
                }

                // remove event from the list (ie shift it to the end)
                ++n_processed;

                // swap with last element that is not yet deactivated
                std::iter_swap(eventIt, events.end() - n_processed);

                // this element's cumulative rate gets initialized with its own rate
                eventIt->cumulative_rate = eventIt->own_rate;
                if (eventIt > events.begin()) {
                    // and then increased by the cumulative rate of its predecessor
                    eventIt->cumulative_rate += (eventIt - 1)->cumulative_rate;
                }
                // now update the cumulative rates of all following elements...
                auto cumsum = (*eventIt).cumulative_rate;
                for (auto _it = eventIt + 1; _it < events.end() - n_processed; ++_it) {
                    cumsum += (*_it).own_rate;
                    (*_it).cumulative_rate = cumsum;
                }

            } else {
                log::critical("this should not happen (event not found by drawn rate)");
                throw std::logic_error("this should not happen (event not found by drawn rate)");
            }
        }

        if(!new_topologies.empty()) {
            for(auto&& top : new_topologies) {
                // if we have a single particle that is not of flavor topology, ignore!
                bool ignore = false;
                if(top.getNParticles() == 1) {
                    auto type = kernel->getCPUKernelStateModel().getParticleData()->entry_at(top.getParticles().front()).type;
                    if(kernel->getKernelContext().particle_types().info_of(type).flavor != readdy::model::Particle::FLAVOR_TOPOLOGY) {
                        ignore = true;
                    }
                }
                if(!ignore) {
                    topologies.push_back(std::make_unique<readdy::model::top::GraphTopology>(std::move(top)));
                    (*topologies.end())->updateReactionRates();
                }
            }
        }

    }

}

}
}
}
}
}