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
 * This header contains the definitions of the threading config, particularly important for the CPU kernel.
 *
 * @file Config.h
 * @brief Config class header
 * @author clonker
 * @date 05.09.16
 */

#pragma once

#include <thread>

#include <readdy/common/common.h>
#include "executor.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(util)
NAMESPACE_BEGIN(thread)

/**
 * Strongly typed enum which contains the different threading modes.
 */
enum class ThreadMode {
    inactive, pool, std_thread, std_async
};

/**
 * Struct that holds the threading configuration, i.e., how many threads should be used when executing code on the
 * CPU kernel.
 */
struct Config {
    /**
     * return type of std::thread::hardware_concurrency()
     */
    using n_threads_type = decltype(std::thread::hardware_concurrency());

    /**
     * constructs a new config (should only be performed by the kernels)
     */
    Config();
     /**
      * constructs a new config with a certain mode
      * @param mode the mode
      */
    Config(ThreadMode mode);

    /**
     * Returns the number of threads. Defaults to:
     *  - hardware_concurrency() if in DEBUG mode
     *  - 4 * hardware_concurrency() otherwise
     * @return the number of threads
     */
    n_threads_type nThreads() const {
        return m_nThreads;
    };

    /**
     * Set the number of threads to be used
     */
    void setNThreads(n_threads_type n) {
        m_nThreads = n;
        update();
    };

    /**
     * Sets the threading mode.
     * @param mode the mode
     */
    void setMode(ThreadMode mode) {
        _mode = mode;
        update();
    };

    /**
     * Yields a pointer to the executor selected by the threading mode.
     * @return a pointer to the configured executor
     */
    const executor_base *const executor() const {
        return _executor.get();
    };

private:
    n_threads_type m_nThreads;
    ThreadMode _mode{ThreadMode::std_thread};
    std::unique_ptr<executor_base> _executor;
    std::unique_ptr<ctpl::thread_pool> pool;

    void update();
};

NAMESPACE_END(thread)
NAMESPACE_END(util)
NAMESPACE_END(readdy)
