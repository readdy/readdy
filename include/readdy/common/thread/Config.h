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

enum class ThreadMode {
    inactive, pool, std_thread, std_async
};

/**
 * Struct that holds the threading configuration, i.e., how many threads should be used when executing code on the
 * CPU or CPU_Dense kernel.
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
     * Returns the number of threads. Defaults to:
     *  - hardware_concurrency() if in DEBUG mode
     *  - 4 * hardware_concurrency() otherwise
     * @return the number of threads
     */
    n_threads_type nThreads() const;

    /**
     * Set the number of threads to be used
     */
    void setNThreads(n_threads_type n);

    void setMode(ThreadMode mode);

    const executor_base *const executor() const;

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
