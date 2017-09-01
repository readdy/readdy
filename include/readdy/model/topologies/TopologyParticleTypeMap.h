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
 * @file TopologyParticleTypeMap.h
 * @brief << brief description >>
 * @author clonker
 * @date 29.08.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <readdy/common/common.h>
#include <readdy/common/hash.h>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(top)

using topology_type_pair = std::tuple<topology_type_type, topology_type_type>;
using topology_particle_type_tuple = std::tuple<particle_type_type, topology_type_type, particle_type_type, topology_type_type>;

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
