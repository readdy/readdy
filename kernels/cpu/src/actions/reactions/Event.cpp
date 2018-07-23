/********************************************************************
 * Copyright © 2018 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * Redistribution and use in source and binary forms, with or       *
 * without modification, are permitted provided that the            *
 * following conditions are met:                                    *
 *  1. Redistributions of source code must retain the above         *
 *     copyright notice, this list of conditions and the            *
 *     following disclaimer.                                        *
 *  2. Redistributions in binary form must reproduce the above      *
 *     copyright notice, this list of conditions and the following  *
 *     disclaimer in the documentation and/or other materials       *
 *     provided with the distribution.                              *
 *  3. Neither the name of the copyright holder nor the names of    *
 *     its contributors may be used to endorse or promote products  *
 *     derived from this software without specific                  *
 *     prior written permission.                                    *
 *                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND           *
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,      *
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF         *
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE         *
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR            *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,     *
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,         *
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER *
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,      *
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)    *
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF      *
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
 ********************************************************************/


/**
 * << detailed description >>
 *
 * @file Event.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 28.10.16
 */

#include <readdy/kernel/cpu/actions/reactions/Event.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace actions {
namespace reactions {

Event::Event(std::uint8_t nEducts, std::uint8_t nProducts, index_type idx1, index_type idx2, scalar reactionRate,
             scalar cumulativeRate, reaction_index_type reactionIdx, ParticleTypeId t1, ParticleTypeId t2)
        : nEducts(nEducts), nProducts(nProducts), idx1(idx1), idx2(idx2), rate(reactionRate),
          cumulativeRate(cumulativeRate), reactionIndex(reactionIdx), t1(t1), t2(t2) {
}

std::ostream &operator<<(std::ostream &os, const Event &evt) {
    os << "Event(" << evt.idx1 << "[type=" << evt.t1 << "]";
    if (evt.nEducts == 2) {
        os << " + " << evt.idx2 << "[type=" << evt.t2 << "]";
    }
    os << ", rate=" << evt.rate << ", cumulativeRate=" << evt.cumulativeRate
       << ", reactionIdx=" << evt.reactionIndex;
    return os;
}

}
}
}
}
}