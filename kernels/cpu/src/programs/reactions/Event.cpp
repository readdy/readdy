/**
 * << detailed description >>
 *
 * @file Event.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 28.10.16
 */

#include <readdy/kernel/cpu/programs/reactions/Event.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace programs {
namespace reactions {

Event::Event(unsigned int nEducts, unsigned int nProducts, index_type idx1, index_type idx2, double reactionRate,
             double cumulativeRate, reaction_index_type reactionIdx, unsigned int t1, unsigned int t2)
        : nEducts(nEducts), nProducts(nProducts), idx1(idx1), idx2(idx2), reactionRate(reactionRate),
          cumulativeRate(cumulativeRate), reactionIdx(reactionIdx), t1(t1), t2(t2) {
}

std::ostream &operator<<(std::ostream &os, const Event &evt) {
    os << "Event(" << evt.idx1 << "[type=" << evt.t1 << "]";
    if (evt.nEducts == 2) {
        os << " + " << evt.idx2 << "[type=" << evt.t2 << "]";
    }
    os << ", rate=" << evt.reactionRate << ", cumulativeRate=" << evt.cumulativeRate
       << ", reactionIdx=" << evt.reactionIdx;
    return os;
}

}
}
}
}
}