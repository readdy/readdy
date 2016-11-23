/**
 * << detailed description >>
 *
 * @file ReactionUtils.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 23.11.16
 */


#include <readdy/kernel/cpu_dense/programs/reactions/ReactionUtils.h>

namespace readdy {
namespace kernel {
namespace cpu_dense {
namespace programs {
namespace reactions {

data_t::update_t handleEventsGillespie(
        Kernel const *const kernel,
        bool filterEventsInAdvance, bool approximateRate,
        std::vector<event_t> &&events) {
    using rdy_particle_t = readdy::model::Particle;

    data_t::update_t newParticles{};

    if(!events.empty()) {
        const auto &ctx = kernel->getKernelContext();
        const auto data = kernel->getKernelStateModel().getParticleData();
        const auto dt = ctx.getTimeStep();
        /**
         * Handle gathered reaction events
         */
        {
            std::size_t nDeactivated = 0;
            const std::size_t nEvents = events.size();
            while (nDeactivated < nEvents) {
                const auto alpha = (*(events.end() - nDeactivated - 1)).cumulativeRate;
                const auto x = readdy::model::rnd::uniform_real(0., alpha);
                const auto eventIt = std::lower_bound(
                        events.begin(), events.end() - nDeactivated, x,
                        [](const event_t &elem1, double elem2) {
                            return elem1.cumulativeRate < elem2;
                        }
                );
                const auto event = *eventIt;
                if (eventIt == events.end() - nDeactivated) {
                    throw std::runtime_error("this should not happen (event not found)");
                }
                if (filterEventsInAdvance || shouldPerformEvent(event.reactionRate, dt, approximateRate)) {
                    /**
                     * Perform reaction
                     */
                    {

                        auto entry1 = event.idx1;
                        if (event.nEducts == 1) {
                            auto reaction = ctx.getOrder1Reactions(event.t1)[event.reactionIdx];
                            performReaction(*data, entry1, entry1, newParticles, reaction);
                        } else {
                            auto reaction = ctx.getOrder2Reactions(event.t1, event.t2)[event.reactionIdx];
                            performReaction(*data, entry1, event.idx2, newParticles, reaction);
                        }
                    }
                    /**
                     * deactivate events whose educts have disappeared (including the just handled one)
                     */
                    {
                        auto _it = events.begin();
                        double cumsum = 0.0;
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
    return std::move(newParticles);
}
}
}
}
}
}