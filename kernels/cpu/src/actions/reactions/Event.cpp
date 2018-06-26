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