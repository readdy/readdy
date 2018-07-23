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
 * This header file contains the definition of joining_future, an RAII implementation of std::future.
 *
 * @file joining_future.h
 * @brief Definition of readdy::util::thread::joining_future
 * @author clonker
 * @date 13.06.17
 * @copyright BSD-3
 */

#pragma once

#include <future>
#include <readdy/common/macros.h>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(util)
NAMESPACE_BEGIN(thread)

template<typename T>
class joining_future {
    std::future<T> _future;

public:
    /**
     * Generates a joining future out of a future
     * @param f the future
     */
    explicit joining_future(std::future<T> &&f) : _future(std::move(f)) {}

    /**
     * move is o.k.
     */
    joining_future(joining_future &&) noexcept = default;

    /**
     * move is o.k.
     */
    joining_future &operator=(joining_future &&) noexcept = default;

    /**
     * no copying
     */
    joining_future(const joining_future &) = delete;

    /**
     * no copying
     */
    joining_future &operator=(const joining_future &) = delete;

    /**
     * wait for the future if it is valid
     */
    ~joining_future() {
        if (_future.valid()) _future.wait();
    }
};

NAMESPACE_END(thread)
NAMESPACE_END(util)
NAMESPACE_END(readdy)
