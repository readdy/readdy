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
 * @file ReactionUtils.h<
 * @brief << brief description >>
 * @author clonker
 * @date 20.10.16
 */

#ifndef READDY_DENSE_KERNEL_REACTIONUTILS_H
#define READDY_DENSE_KERNEL_REACTIONUTILS_H

#include <cmath>
#include <readdy/model/RandomProvider.h>
#include <readdy/common/logging.h>
#include "Event.h"
#include "readdy/kernel/cpu_dense/CPUDKernel.h"

namespace readdy {
namespace kernel {
namespace cpu_dense {
namespace programs {
namespace reactions {

using kernel_t = readdy::kernel::cpu_dense::CPUDKernel;
using vec_t = readdy::model::Vec3;
using data_t = readdy::kernel::cpu_dense::model::CPUDParticleData;
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

data_t::update_t handleEventsGillespie(
        CPUDKernel const *const kernel,
        bool filterEventsInAdvance, bool approximateRate,
        std::vector<event_t> &&events);

template<typename ParticleIndexCollection>
void gatherEvents(CPUDKernel const *const kernel, const ParticleIndexCollection &particles, const nl_t &nl,
                  const data_t &data, double &alpha, std::vector<event_t> &events,
                  const readdy::model::KernelContext::dist_squared_fun& d2) {
    for (const auto index : particles) {
        auto& entry = data.entry_at(index);
        // this being false should really not happen, though
        if (!entry.deactivated) {
            // order 1
            {
                const auto &reactions = kernel->getKernelContext().getOrder1Reactions(entry.type);
                for (auto it = reactions.begin(); it != reactions.end(); ++it) {
                    const auto rate = (*it)->getRate();
                    if (rate > 0) {
                        alpha += rate;
                        events.push_back(
                                {1, (*it)->getNProducts(), index, 0, rate, alpha,
                                 static_cast<event_t::reaction_index_type>(it - reactions.begin()),
                                 entry.type, 0});
                    }
                }
            }
            // order 2
            for (const auto& idx_neighbor : nl->find_neighbors(index)) {
                if (index > idx_neighbor.idx) continue;
                const auto& neighbor = data.entry_at(idx_neighbor.idx);
                const auto &reactions = kernel->getKernelContext().getOrder2Reactions(entry.type, neighbor.type);
                if (!reactions.empty()) {
                    const auto distSquared = d2(neighbor.position(), entry.position());
                    for (auto it = reactions.begin(); it < reactions.end(); ++it) {
                        const auto &react = *it;
                        const auto rate = react->getRate();
                        if (rate > 0 && distSquared < react->getEductDistanceSquared()) {
                            alpha += rate;
                            events.push_back({2, react->getNProducts(), index, idx_neighbor.idx,
                                              rate, alpha,
                                              static_cast<event_t::reaction_index_type>(it - reactions.begin()),
                                              entry.type, neighbor.type});
                        }
                    }
                }
            }
        }
    }

}

template<bool checkDeactivated=false, typename Reaction>
void performReaction(data_t& data, data_t::index_t idx1, data_t::index_t idx2, data_t::update_t& newEntries, Reaction* reaction) {
    if(checkDeactivated) {
        if(data.entry_at(idx1).deactivated || data.entry_at(idx2).deactivated) return;
    }
    switch(reaction->getType()) {
        case reaction_type::Decay: {
            data.deactivate(idx1);
            break;
        }
        case reaction_type::Conversion: {
            data.entry_at(idx1).type = reaction->getProducts()[0];
            break;
        }
        case reaction_type::Enzymatic: {
            auto& entry1 = data.entry_at(idx1);
            if (entry1.type == reaction->getEducts()[1]) {
                // p1 is the catalyst
                data.entry_at(idx2).type = reaction->getProducts()[0];
            } else {
                // p2 is the catalyst
                entry1.type = reaction->getProducts()[0];
            }
            break;
        }
        case reaction_type::Fission: {
            auto& entry1 = data.entry_at(idx1);
            auto n3 = readdy::model::rnd::normal3(0, 1);
            n3 /= std::sqrt(n3 * n3);

            readdy::model::Particle p (entry1.position() - reaction->getWeight2() * reaction->getProductDistance() * n3, reaction->getProducts()[1]);
            newEntries.push_back({p});

            entry1.type = reaction->getProducts()[0];
            data.displace(entry1, reaction->getWeight1() * reaction->getProductDistance() * n3);
            break;
        }
        case reaction_type::Fusion: {
            auto& entry1 = data.entry_at(idx1);
            auto& entry2 = data.entry_at(idx2);
            const auto& e1Pos = entry1.position();
            const auto& e2Pos = entry2.position();
            if (reaction->getEducts()[0] == entry1.type) {
                data.displace(entry1, reaction->getWeight1() * (e2Pos - e1Pos));
            } else {
                data.displace(entry1, reaction->getWeight1() * (e1Pos - e2Pos));
            }
            entry1.type = reaction->getProducts()[0];
            data.deactivate(entry2);
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
