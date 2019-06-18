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
 * @file TopologyParticleTypeMap.h
 * @brief << brief description >>
 * @author clonker
 * @date 29.08.17
 * @copyright BSD-3
 */

#pragma once

#include <readdy/common/common.h>
#include <readdy/common/hash.h>

namespace readdy::model::top {

using topology_type_pair = std::tuple<TopologyTypeId, TopologyTypeId>;
using topology_particle_type_tuple = std::tuple<ParticleTypeId, TopologyTypeId, ParticleTypeId, TopologyTypeId>;

namespace detail {

class TopologyParticleTypeHasher {
public:
    /**
     * This hasher needs to respect that in the tuple (x1, t1, x2, t2),
     * swapping (x1, t1) with (x2, t2) should result in the same hash.
     * Thus first compute hashes of both sub-tuples and compare them accordingly, which already gives a unique hash.
     * A weaker/faster version e.g. could only hash (t1, t2) and neglect x1 and x2, leaving the rest to the comparator.
     */
    std::size_t operator()(const topology_particle_type_tuple &tup) const {
        std::size_t seed1{0};
        util::hash::combine(seed1, std::get<0>(tup));
        util::hash::combine(seed1, std::get<1>(tup));
        std::size_t seed2{0};
        util::hash::combine(seed2, std::get<2>(tup));
        util::hash::combine(seed2, std::get<3>(tup));

        std::size_t seed{0};
        if (seed1 > seed2) {
            util::hash::combine(seed, seed1);
            util::hash::combine(seed, seed2);
        } else {
            util::hash::combine(seed, seed2);
            util::hash::combine(seed, seed1);
        }
        return seed;
    }
};

class TopologyParticleTypeEq {
public:
    bool operator()(const topology_particle_type_tuple &t1, const topology_particle_type_tuple &t2) const {
        return t1 == t2
               || (std::get<0>(t1) == std::get<2>(t2) && std::get<1>(t1) == std::get<3>(t2)
                   && std::get<2>(t1) == std::get<0>(t2) && std::get<3>(t1) == std::get<1>(t2));
    }
};

}

template<typename T>
using topology_particle_type_tuple_umap = std::unordered_map<topology_particle_type_tuple, T,
                                                             detail::TopologyParticleTypeHasher,
                                                             detail::TopologyParticleTypeEq>;

}
