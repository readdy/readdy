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

#pragma once
#include <spdlog/fmt/ostr.h>
#include <readdy/model/Particle.h>
#include <readdy/model/reactions/Reaction.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace actions {
namespace reactions {
struct Event {
    using index_type = std::size_t;
    using reaction_index_type = std::size_t;
    std::uint8_t nEducts;
    std::uint8_t nProducts;
    index_type idx1, idx2;
    reaction_index_type reactionIndex;
    readdy::model::Particle::type_type t1, t2;
    readdy::scalar rate;
    readdy::scalar cumulativeRate;

    Event(std::uint8_t nEducts, std::uint8_t nProducts, index_type idx1, index_type idx2, readdy::scalar reactionRate,
          readdy::scalar cumulativeRate, reaction_index_type reactionIdx, ParticleTypeId t1, ParticleTypeId t2);

    friend std::ostream &operator<<(std::ostream &, const Event &);

};
}
}
}
}
}
