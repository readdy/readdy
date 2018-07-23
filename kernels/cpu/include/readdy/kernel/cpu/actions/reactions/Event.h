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
