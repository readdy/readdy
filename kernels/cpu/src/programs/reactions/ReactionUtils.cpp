/**
 * << detailed description >>
 *
 * @file ReactionUtils.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 20.10.16
 */

#include <readdy/kernel/cpu/programs/reactions/ReactionUtils.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace programs {
namespace reactions {

data_t::entries_t handleEventsGillespie(
        CPUKernel const *const kernel,
        bool filterEventsInAdvance, bool approximateRate,
        std::vector<readdy::kernel::singlecpu::programs::reactions::ReactionEvent> &&events) {
    using event_t = readdy::kernel::singlecpu::programs::reactions::ReactionEvent;
    using rdy_particle_t = readdy::model::Particle;

    data_t::entries_t newParticles{};

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

                    auto& entry1 = data->entries[event.idx1];
                    if (event.nEducts == 1) {
                        auto reaction = ctx.getOrder1Reactions(event.t1)[event.reactionIdx];
                        performReaction(*data, entry1, entry1, event.idx1, event.idx2, newParticles, reaction);
                    } else {
                        auto reaction = ctx.getOrder2Reactions(event.t1, event.t2)[event.reactionIdx];
                        auto& entry2 = data->entries[event.idx2];
                        performReaction(*data, entry1, entry2, event.idx1, event.idx2, newParticles, reaction);
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
    return newParticles;
}
}
}
}
}
}