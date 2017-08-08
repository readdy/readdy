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


CPUEvaluateTopologyReactions::CPUEvaluateTopologyReactions(CPUKernel *const kernel, scalar timeStep)
        : EvaluateTopologyReactions(timeStep), kernel(kernel) {}

template<bool approximated>
bool performReactionEvent(const scalar rate, const scalar timeStep) {
    if (approximated) {
        return readdy::model::rnd::uniform_real() < rate * timeStep;
    } else {
        return readdy::model::rnd::uniform_real() < 1 - std::exp(-rate * timeStep);
    }
}

bool shouldPerformEvent(const scalar rate, const scalar timestep, bool approximated) {
    return approximated ? performReactionEvent<true>(rate, timestep) : performReactionEvent<false>(rate, timestep);
}

void CPUEvaluateTopologyReactions::perform() {
    using rate_t = readdy::model::top::GraphTopology::rate_t;
    auto &topologies = kernel->getCPUKernelStateModel().topologies();

    if(!topologies.empty()) {
        std::stringstream ss;
        for(const auto& top : topologies) {
            if(!top->isDeactivated()) {
                const auto * address = static_cast<const void*>(top.get());
                std::stringstream ss2;
                ss2 << address;
                std::string name = ss2.str();
                ss << ", " << name << "(" << top->getNParticles() <<")";
            }
        }
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

                        TREvent event{};
                        event.own_rate = std::get<1>(reaction);
                        event.cumulative_rate = event.own_rate + current_cumulative_rate;
                        current_cumulative_rate = event.cumulative_rate;
                        event.topology_idx = topology_idx;
                        event.reaction_idx = reaction_idx;

                        events.push_back(event);

                        ++reaction_idx;
                    }
                }
                ++topology_idx;
            }
        }

        if (!events.empty()) {

            auto end = events.end();

            std::vector<readdy::model::top::GraphTopology> new_topologies;
            while (end != events.begin()) {
                const auto cumulative_rate = (end - 1)->cumulative_rate;

                const auto x = readdy::model::rnd::uniform_real(static_cast<scalar>(0.), cumulative_rate);

                const auto eventIt = std::lower_bound(
                        events.begin(), end, x, [](const TREvent &elem1, const rate_t elem2) {
                            return elem1.cumulative_rate < elem2;
                        }
                );

                if (eventIt != events.end()) {
                    const auto &event = *eventIt;

                    if (shouldPerformEvent(event.own_rate, timeStep, true)) {
                        log::trace("picked event {} / {} with rate {}", std::distance(events.begin(), eventIt)+1, events.size(), eventIt->own_rate);
                        // perform the event!
                        auto &topology = topologies.at(event.topology_idx);
                        if (topology->isDeactivated()) {
                            log::critical("deactivated topology with idx {}", event.topology_idx);
                            for (auto it = events.begin(); it != end; ++it) {
                                log::warn(" -> event {}, {}, {}, {}", it->topology_idx, it->reaction_idx, it->own_rate,
                                          it->cumulative_rate);
                            }
                            for (auto it = end; it != events.end(); ++it) {
                                log::warn(" -> deactivated event {}, {}, {}, {}", it->topology_idx, it->reaction_idx,
                                          it->own_rate, it->cumulative_rate);
                            }
                        }
                        assert(!topology->isDeactivated());
                        auto &reaction = topology->registeredReactions().at(event.reaction_idx);
                        auto result = std::get<0>(reaction).execute(*topology, kernel);
                        if (!result.empty()) {
                            // we had a topology fission, so we need to actually remove the current topology from the
                            // data structure
                            topologies.erase(topologies.begin() + event.topology_idx);
                            //log::error("erased topology with index {}", event.topology_idx);
                            assert(topology->isDeactivated());
                            std::move(result.begin(), result.end(), std::back_inserter(new_topologies));
                        } else {
                            if (topology->isNormalParticle(*kernel)) {
                                topologies.erase(topologies.begin() + event.topology_idx);
                                //log::error("erased topology with index {}", event.topology_idx);
                                assert(topology->isDeactivated());
                            }
                        }

                        // loop over all other events, swap events considering this particular topology to the end
                        rate_t cumsum = 0;
                        //log::warn("------");
                        for (auto it = events.begin();;) {
                            if (it->topology_idx == event.topology_idx) {
                                --end;
                                //log::warn("swapping event with topology idx {}", it->topology_idx);
                                std::iter_swap(it, end);
                            }
                            cumsum += it->own_rate;
                            it->cumulative_rate = cumsum;
                            if (it >= end) {
                                break;
                            }
                            if (it->topology_idx != event.topology_idx) {
                                ++it;
                            }
                        }

                    } else {

                        // remove event from the list (ie shift it to the end)
                        --end;

                        // swap with last element that is not yet deactivated
                        std::iter_swap(eventIt, end);

                        // this element's cumulative rate gets initialized with its own rate
                        eventIt->cumulative_rate = eventIt->own_rate;
                        if (eventIt > events.begin()) {
                            // and then increased by the cumulative rate of its predecessor
                            eventIt->cumulative_rate += (eventIt - 1)->cumulative_rate;
                        }
                        // now update the cumulative rates of all following elements...
                        auto cumsum = (*eventIt).cumulative_rate;
                        for (auto _it = eventIt + 1; _it < end; ++_it) {
                            cumsum += (*_it).own_rate;
                            (*_it).cumulative_rate = cumsum;
                        }

                    }

                } else {
                    log::critical("this should not happen (event not found by drawn rate)");
                    throw std::logic_error("this should not happen (event not found by drawn rate)");
                }
            }

            if (!new_topologies.empty()) {
                for (auto &&top : new_topologies) {
                    // if we have a single particle that is not of flavor topology, ignore!
                    if (!top.isNormalParticle(*kernel)) {
                        auto new_top = std::make_unique<readdy::model::top::GraphTopology>(std::move(top));
                        auto it = topologies.push_back(std::move(new_top));
                        (*it)->updateReactionRates();
                        (*it)->configure();
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
}