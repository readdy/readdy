/**
 * << detailed description >>
 *
 * @file ReactionUtils.h
 * @brief << brief description >>
 * @author clonker
 * @date 20.10.16
 */

#ifndef READDY_CPUKERNEL_REACTIONUTILS_H
#define READDY_CPUKERNEL_REACTIONUTILS_H

#include <cmath>
#include <readdy/model/RandomProvider.h>
#include <readdy/kernel/singlecpu/programs/SingleCPUReactionImpls.h>
#include <readdy/kernel/cpu/CPUKernel.h>
#include <readdy/common/logging.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace programs {
namespace reactions {

using kernel_t = readdy::kernel::cpu::CPUKernel;
using vec_t = readdy::model::Vec3;
using data_t = decltype(std::declval<kernel_t>().getKernelStateModel().getParticleData());
using nl_t = const decltype(std::declval<kernel_t>().getKernelStateModel().getNeighborList());
using ctx_t = std::remove_const<decltype(std::declval<kernel_t>().getKernelContext())>::type;
using event_t = readdy::kernel::singlecpu::programs::reactions::ReactionEvent;
using event_index_t = event_t::index_type;
using particle_t = readdy::model::Particle;

template<bool approximated>
bool performReactionEvent(const double rate, const double timeStep) {
    if (approximated) {
        return readdy::model::rnd::uniform_real() < rate * timeStep;
    } else {
        return readdy::model::rnd::uniform_real() < 1 - std::exp(-rate * timeStep);
    }
}


inline bool shouldPerformEvent(const double rate, const double timestep, bool approximated) {
    return approximated ? performReactionEvent<true>(rate, timestep) : performReactionEvent<false>(rate, timestep);
}

std::vector<readdy::model::Particle> handleEventsGillespie(
        CPUKernel const *const kernel,
        bool filterEventsInAdvance, bool approximateRate,
        std::vector<readdy::kernel::singlecpu::programs::reactions::ReactionEvent> &&events);

template<typename ParticleIndexCollection>
void gatherEvents(CPUKernel const *const kernel, const ParticleIndexCollection &particles, const nl_t &nl,
                  const data_t &data, double &alpha, std::vector<event_t> &events) {
    for (const auto idx : particles) {
        // this being false should really not happen, though
        if (!*(data->begin_deactivated() + idx)) {
            const auto particleType = *(data->begin_types() + idx);
            // order 1
            {
                const auto &reactions = kernel->getKernelContext().getOrder1Reactions(particleType);
                for (auto it = reactions.begin(); it != reactions.end(); ++it) {
                    const auto rate = (*it)->getRate();
                    if (rate > 0) {
                        alpha += rate;
                        events.push_back(
                                {1, (*it)->getNProducts(), idx, 0, rate, alpha,
                                 static_cast<event_index_t>(it - reactions.begin()),
                                 particleType, 0});
                    }
                }
            }
            // order 2
            {
                auto nl_it = nl->pairs->find(idx);
                if (nl_it != nl->pairs->end()) {
                    for (const auto &idx_neighbor : nl_it->second) {
                        if (idx > idx_neighbor.idx) continue;
                        const auto neighborType = *(data->begin_types() + idx_neighbor.idx);
                        const auto &reactions = kernel->getKernelContext().getOrder2Reactions(
                                particleType, neighborType
                        );
                        if (!reactions.empty()) {
                            const auto distSquared = idx_neighbor.d2;
                            for (auto it = reactions.begin(); it < reactions.end(); ++it) {
                                const auto &react = *it;
                                const auto rate = react->getRate();
                                if (rate > 0 && distSquared < react->getEductDistanceSquared()) {
                                    alpha += rate;
                                    events.push_back({2, react->getNProducts(), idx, idx_neighbor.idx,
                                                      rate, alpha,
                                                      static_cast<event_index_t>(it - reactions.begin()),
                                                      particleType, neighborType});
                                }
                            }
                        }
                    }
                }
            }
        } else {
            log::console()->error("The particles list which was given to gather events contained a particle that "
                                          "was already deactivated. This should not happen!");
        }
    }

}

}
}
}
}
}
#endif //READDY_CPUKERNEL_REACTIONUTILS_H
