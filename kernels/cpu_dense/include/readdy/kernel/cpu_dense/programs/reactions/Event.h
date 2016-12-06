/**
 * << detailed description >>
 *
 * @file ReactionEvent.h
 * @brief << brief description >>
 * @author clonker
 * @date 22.11.16
 */

#ifndef READDY_MAIN_REACTIONEVENT_H
#define READDY_MAIN_REACTIONEVENT_H

#include <spdlog/fmt/ostr.h>
#include <readdy/kernel/cpu_dense/model/CPUDParticleData.h>

namespace readdy {
namespace kernel {
namespace cpu_dense {
namespace programs {
namespace reactions {
struct Event {
    using index_type = model::CPUDParticleData::index_t;
    using reaction_index_type = std::size_t;
    unsigned int nEducts;
    unsigned int nProducts;
    index_type idx1, idx2;
    reaction_index_type reactionIdx;
    unsigned int t1, t2;
    double reactionRate;
    double cumulativeRate;

    Event(unsigned int nEducts, unsigned int nProducts, index_type idx1, index_type idx2, double reactionRate,
          double cumulativeRate, reaction_index_type reactionIdx, unsigned int t1, unsigned int t2);

    friend std::ostream &operator<<(std::ostream &, const Event &);

};
}
}
}
}
}
#endif //READDY_MAIN_REACTIONEVENT_H
