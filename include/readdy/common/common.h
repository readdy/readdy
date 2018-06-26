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
 * @file common.h
 * @brief << brief description >>
 * @author clonker
 * @date 07.03.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include "logging.h"
#include "ReaDDyVec3.h"

NAMESPACE_BEGIN(h5rd)
class File;
class Group;
NAMESPACE_END(h5rd)

NAMESPACE_BEGIN(readdy)

constexpr inline std::size_t operator "" _z ( unsigned long long n ) { return n; }

using scalar = double;
using stride_type = std::uint32_t;
using Vec3 = _internal::ReaDDyVec3<scalar>;
using Matrix33 = _internal::ReaDDyMatrix33<scalar>;
using time_step_type = unsigned long;
using ParticleTypeId = unsigned short;
// signed short as sometimes it is needed to have a "blank" topology type (corresponding to -1)
using TopologyTypeId = short;
constexpr TopologyTypeId EmptyTopologyId = static_cast<TopologyTypeId>(-1);

constexpr bool single_precision = std::is_same<scalar, float>::value;
constexpr bool double_precision = std::is_same<scalar, double>::value;

NAMESPACE_BEGIN(c_)
constexpr scalar zero = static_cast<scalar>(0.0);
constexpr scalar one = static_cast<scalar>(1.0);
constexpr scalar two = static_cast<scalar>(2.0);
constexpr scalar three = static_cast<scalar>(3.0);
constexpr scalar four = static_cast<scalar>(4.0);
constexpr scalar five = static_cast<scalar>(5.0);
constexpr scalar half = static_cast<scalar>(.5);
NAMESPACE_END(c_)

using File = h5rd::File;

inline auto readdy_default_n_threads() -> decltype(std::thread::hardware_concurrency()) {
    using return_type = decltype(std::thread::hardware_concurrency());
#ifdef READDY_DEBUG
    return_type m_nThreads = std::max(std::thread::hardware_concurrency(), 1u);
#else
    // magic number 4 to enable some load balancing
    return_type m_nThreads = std::max(4*std::thread::hardware_concurrency(), 1u);
#endif
    const char *env = std::getenv("READDY_N_CORES");
    if (env != nullptr) {
        m_nThreads = static_cast<decltype(m_nThreads)>(std::stol(env));
        log::debug("Using {} threads (set by environment variable READDY_N_CORES)", m_nThreads);
    }
    return m_nThreads;
}

NAMESPACE_END(readdy)
