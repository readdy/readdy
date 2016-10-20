/**
 * << detailed description >>
 *
 * @file ReactionUtils.h
 * @brief << brief description >>
 * @author clonker
 * @date 20.10.16
 */

#ifndef READDY_MAIN_REACTIONUTILS_H
#define READDY_MAIN_REACTIONUTILS_H

#include <cmath>
#include <readdy/model/RandomProvider.h>
#include <readdy/kernel/singlecpu/programs/SingleCPUReactionImpls.h>
#include <readdy/kernel/cpu/CPUKernel.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace programs {
namespace reactions {

template<bool approximated>
bool performReactionEvent(const double rate, const double timeStep) {
    if (approximated) {
        return readdy::model::rnd::uniform() < rate * timeStep;
    } else {
        return readdy::model::rnd::uniform() < 1 - std::exp(-rate * timeStep);
    }
}


inline bool performEvent(const double rate, const double timestep, bool approximated) {
    return approximated ? performReactionEvent<true>(rate, timestep) : performReactionEvent<false>(rate, timestep);
}

std::vector<readdy::model::Particle> handleEventsGillespie(
        CPUKernel const *const kernel,
        bool filterEventsInAdvance, bool approximateRate,
        std::vector<readdy::kernel::singlecpu::programs::reactions::ReactionEvent> &&events);

}
}
}
}
}
#endif //READDY_MAIN_REACTIONUTILS_H
