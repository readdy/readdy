/**
 * << detailed description >>
 *
 * @file Gillespie.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 20.10.16
 */
#include <readdy/kernel/cpu/programs/reactions/Gillespie.h>


namespace readdy {
namespace kernel {
namespace cpu {
namespace programs {
namespace reactions {

Gillespie::Gillespie(const CPUKernel *const kernel) : kernel(kernel) {}

std::vector<singlecpu::programs::reactions::ReactionEvent> Gillespie::gatherEvents(double &alpha) {
    using index_t = singlecpu::programs::reactions::ReactionEvent::index_type;
    std::vector<singlecpu::programs::reactions::ReactionEvent> events;
    const auto &ctx = kernel->getKernelContext();
    auto data = kernel->getKernelStateModel().getParticleData();
    /**
    * Reactions with one educt
    */
    {
        auto it_type = data->begin_types();
        auto it_deactivated = data->begin_deactivated();
        const auto end = data->end_types();
        while (it_type != end) {
            if (!*it_deactivated) {
                const auto &reactions = ctx.getOrder1Reactions(*it_type);
                for (auto it = reactions.begin(); it != reactions.end(); ++it) {
                    const auto rate = (*it)->getRate();
                    if (rate > 0) {
                        alpha += rate;
                        events.push_back({1, (*it)->getNProducts(),
                                          (index_t) (it_type - data->begin_types()), 0, rate, alpha,
                                          (index_t) (it - reactions.begin()), *it_type, 0});
                    }
                }
            }
            ++it_type;
            ++it_deactivated;
        }
    }

    /**
     * Reactions with two educts
     */
    {
        const auto *neighborList = kernel->getKernelStateModel().getNeighborList();
        auto typesBegin = data->begin_types();
        for (auto &&nl_it = neighborList->pairs->begin(); nl_it != neighborList->pairs->end(); ++nl_it) {
            const index_t idx1 = nl_it->first;
            if (*(data->begin_deactivated() + idx1)) continue;
            auto neighbors = nl_it->second;
            for (const auto &neighbor : neighbors) {
                if (idx1 > neighbor.idx) continue;
                if (*(data->begin_deactivated() + neighbor.idx)) continue;
                const auto &reactions = ctx.getOrder2Reactions(
                        *(data->begin_types() + idx1), *(data->begin_types() + neighbor.idx)
                );

                const auto distSquared = neighbor.d2;

                for (auto it = reactions.begin(); it != reactions.end(); ++it) {
                    // if close enough
                    const auto reaction = *it;
                    if (distSquared < reaction->getEductDistanceSquared()) {
                        const auto rate = reaction->getRate();
                        if (rate > 0) {
                            alpha += rate;
                            events.push_back(
                                    {2, reaction->getNProducts(), idx1, neighbor.idx, rate, alpha,
                                     (index_t) (it - reactions.begin()),
                                     *(typesBegin + idx1), *(typesBegin + neighbor.idx)});
                        }
                    }
                }
            }
        }
    }
    return events;
}

}
}
}
}
}