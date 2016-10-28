/**
 * << detailed description >>
 *
 * @file ReactionUtils.h<
 * @brief << brief description >>
 * @author clonker
 * @date 20.10.16
 */

#ifndef READDY_CPUKERNEL_REACTIONUTILS_H
#define READDY_CPUKERNEL_REACTIONUTILS_H

#include <cmath>
#include <readdy/model/RandomProvider.h>
#include <readdy/kernel/cpu/CPUKernel.h>
#include <readdy/common/logging.h>
#include "Event.h"

namespace readdy {
namespace kernel {
namespace cpu {
namespace programs {
namespace reactions {

using kernel_t = readdy::kernel::cpu::CPUKernel;
using vec_t = readdy::model::Vec3;
using data_t = readdy::kernel::cpu::model::ParticleData;
using reaction_type = readdy::model::reactions::Reaction<1>::ReactionType;
using nl_t = const decltype(std::declval<kernel_t>().getKernelStateModel().getNeighborList());
using ctx_t = std::remove_const<decltype(std::declval<kernel_t>().getKernelContext())>::type;
using event_t = Event;
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

std::pair<data_t::entries_t, std::vector<data_t::Entry*>> handleEventsGillespie(
        CPUKernel const *const kernel,
        bool filterEventsInAdvance, bool approximateRate,
        std::vector<event_t> &&events);

template<typename ParticleIndexCollection>
void gatherEvents(CPUKernel const *const kernel, const ParticleIndexCollection &particles, const nl_t &nl,
                  const data_t &data, double &alpha, std::vector<event_t> &events) {
    for (const auto idx : particles) {
        // this being false should really not happen, though
        if (!idx->is_deactivated()) {
            // order 1
            {
                const auto &reactions = kernel->getKernelContext().getOrder1Reactions(idx->type);
                for (auto it = reactions.begin(); it != reactions.end(); ++it) {
                    const auto rate = (*it)->getRate();
                    if (rate > 0) {
                        alpha += rate;
                        events.push_back(
                                {1, (*it)->getNProducts(), idx, 0, rate, alpha,
                                 static_cast<event_t::reaction_index_type>(it - reactions.begin()),
                                 idx->type, 0});
                    }
                }
            }
            // order 2
            {
                auto nl_it = nl->pairs.find(idx);
                if (nl_it != nl->pairs.end()) {
                    for (const auto idx_neighbor : nl_it->second) {
                        if (idx > idx_neighbor.idx) continue;
                        const auto neighbor = idx_neighbor.idx;
                        const auto &reactions = kernel->getKernelContext().getOrder2Reactions(idx->type, neighbor->type);
                        if (!reactions.empty()) {
                            const auto distSquared = idx_neighbor.d2;
                            for (auto it = reactions.begin(); it < reactions.end(); ++it) {
                                const auto &react = *it;
                                const auto rate = react->getRate();
                                if (rate > 0 && distSquared < react->getEductDistanceSquared()) {
                                    alpha += rate;
                                    events.push_back({2, react->getNProducts(), idx, idx_neighbor.idx,
                                                      rate, alpha,
                                                      static_cast<event_t::reaction_index_type>(it - reactions.begin()),
                                                      idx->type, neighbor->type});
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

template<typename Reaction>
void performReaction(data_t& data, data_t::Entry* e1, data_t::Entry* e2, data_t::entries_t& newEntries,
                     std::vector<data_t::Entry*>& decayedEntries, Reaction* reaction) {
    switch(reaction->getType()) {
        case reaction_type::Decay: {
            decayedEntries.push_back(e1);
            //data.removeEntry(e1);
            break;
        }
        case reaction_type::Conversion: {
            e1->type = reaction->getProducts()[0];
            break;
        }
        case reaction_type::Enzymatic: {
            if (e1->type == reaction->getEducts()[1]) {
                // p1 is the catalyst
                e2->type = reaction->getProducts()[0];
            } else {
                // p2 is the catalyst
                e1->type = reaction->getProducts()[0];
            }
            break;
        }
        case reaction_type::Fission: {
            auto n3 = readdy::model::rnd::normal3(0, 1);
            n3 /= sqrt(n3 * n3);
            e1->type = reaction->getProducts()[0];
            e1->pos = e1->pos + reaction->getWeight1() * reaction->getProductDistance() * n3;

            readdy::model::Particle p (e1->pos - reaction->getWeight2() * reaction->getProductDistance() * n3, reaction->getProducts()[1]);
            newEntries.push_back({p});
            break;
        }
        case reaction_type::Fusion: {
            e1->type = reaction->getProducts()[0];
            if (reaction->getEducts()[0] == e1->type) {
                e1->pos = e1->pos + reaction->getWeight1() * (e2->pos - e1->pos);
            } else {
                e1->pos = e2->pos + reaction->getWeight1() * (e1->pos - e2->pos);
            }
            decayedEntries.push_back(e2);
            // data.removeEntry(e2);
            break;
        }
    }
}

}
}
}
}
}
#endif //READDY_CPUKERNEL_REACTIONUTILS_H
