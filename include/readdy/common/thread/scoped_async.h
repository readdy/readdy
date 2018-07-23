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
 * The definition of scoped_async. Behaves the same as scoped_thread, just with std::async as backing call.
 *
 * @file scoped_async.h
 * @brief Contains the definition of `readdy::util::thread::scoped_async`
 * @author clonker
 * @date 28.03.17
 * @copyright BSD-3
 */

#pragma once

#include <future>
#include <readdy/common/macros.h>
#include <readdy/common/logging.h>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(util)
NAMESPACE_BEGIN(thread)

class scoped_async {
    std::future<void> async_;
public:
    /**
     * Creates a new scoped_async object that will execute a provided function with the respective arguments
     * asynchronously.
     *
     * @tparam Function the function type
     * @tparam Args the argument types
     * @param fun the function instance
     * @param args the arguments
     */
    template<typename Function, typename... Args>
    explicit scoped_async(Function &&fun, Args &&... args)
            : async_(std::async(std::launch::async, std::forward<Function>(fun), std::forward<Args>(args)...)) {}

    /**
     * wait for the task to finish if valid
     */
    ~scoped_async() {
        if(async_.valid()) async_.wait();
    }

    /**
     * no copying
     */
    scoped_async(const scoped_async &) = delete;

    /**
     * no copy assign
     */
    scoped_async &operator=(const scoped_async &) = delete;

    /**
     * move is permitted
     */
    scoped_async(scoped_async &&) = default;

    /**
     * move assign is permitted
     */
    scoped_async &operator=(scoped_async &&) = default;
};

NAMESPACE_END(thread)
NAMESPACE_END(util)
NAMESPACE_END(readdy)
