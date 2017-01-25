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
 * @file ReactionEvent.h
 * @brief << brief description >>
 * @author clonker
 * @date 28.10.16
 */

#ifndef READDY_MAIN_REACTIONEVENT_H
#define READDY_MAIN_REACTIONEVENT_H

#include <spdlog/fmt/ostr.h>
#include <readdy/kernel/cpu/model/CPUParticleData.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace actions {
namespace reactions {
struct Event {
    using index_type = model::CPUParticleData::index_t;
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
