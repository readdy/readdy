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
 * @file common.h
 * @brief << brief description >>
 * @author clonker
 * @date 07.03.17
 * @copyright BSD-3
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
