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

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(top)

using topology_type_pair = std::tuple<TopologyTypeId, TopologyTypeId>;
using topology_particle_type_tuple = std::tuple<ParticleTypeId, TopologyTypeId, ParticleTypeId, TopologyTypeId>;

NAMESPACE_BEGIN(detail)

class TopologyParticleTypeHasher {
public:
    std::size_t operator()(const topology_particle_type_tuple& tup) const {
        std::size_t seed {0};
        if(std::get<1>(tup) > std::get<3>(tup)) {
            util::hash::combine(seed, std::get<0>(tup));
            util::hash::combine(seed, std::get<1>(tup));
            util::hash::combine(seed, std::get<2>(tup));
            util::hash::combine(seed, std::get<3>(tup));
        } else {
            util::hash::combine(seed, std::get<2>(tup));
            util::hash::combine(seed, std::get<3>(tup));
            util::hash::combine(seed, std::get<0>(tup));
            util::hash::combine(seed, std::get<1>(tup));
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

NAMESPACE_END(detail)

template<typename T>
using topology_particle_type_tuple_umap = std::unordered_map<topology_particle_type_tuple, T,
                                                             detail::TopologyParticleTypeHasher,
                                                             detail::TopologyParticleTypeEq>;

NAMESPACE_END(top)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
