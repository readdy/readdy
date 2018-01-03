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
