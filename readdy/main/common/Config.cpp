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
 * @file Config.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 05.09.16
 */

#include <thread>
#include <algorithm>

#include <readdy/common/logging.h>

#include "readdy/common/thread/Config.h"

#if READDY_OSX
#include <cstdlib>
#endif

namespace readdy {
namespace util {
namespace thread {

Config::Config() {
#ifdef READDY_DEBUG
    m_nThreads = std::max(std::thread::hardware_concurrency(), 1u);
#else
    // magic number 4 to enable some load balancing
    m_nThreads = std::max(4*std::thread::hardware_concurrency(), 1u);
#endif

    const char *env = std::getenv("READDY_N_CORES");
    if (env) {
        m_nThreads = static_cast<n_threads_t>(std::stol(env));
        log::debug("Using {} threads (set by environment variable READDY_N_CORES)", m_nThreads);
    }
    update();
}

Config::n_threads_t Config::nThreads() const {
    return m_nThreads;
}

void Config::setNThreads(const Config::n_threads_t n) {
    m_nThreads = n;
    update();
}

void Config::setMode(ThreadMode mode) {
    _mode = mode;
    update();
}

void Config::update() {
    using pool_executor = readdy::util::thread::executor<readdy::util::thread::executor_type::pool>;
    using thread_executor = readdy::util::thread::executor<readdy::util::thread::executor_type::std_thread>;
    using async_executor = readdy::util::thread::executor<readdy::util::thread::executor_type::std_async>;

    switch (_mode) {
        case ThreadMode::pool: {
            if (pool) {
                pool->resize(m_nThreads);
            } else {
                pool = std::make_unique<ctpl::thread_pool>(m_nThreads);
            }
            _executor = std::unique_ptr<executor_base>(new pool_executor(pool.get()));
            break;
        }
        case ThreadMode::std_thread: {
            if(pool) {
                pool->stop(true);
                pool.reset(nullptr);
            }
            _executor = std::unique_ptr<executor_base>(new thread_executor(nullptr));
            break;
        }
        case ThreadMode::std_async: {
            if(pool) {
                pool->stop(true);
                pool.reset(nullptr);
            }
            _executor = std::unique_ptr<executor_base>(new async_executor(nullptr));
            break;
        }
        case ThreadMode::inactive: {
            if(pool) {
                pool->stop(true);
                pool.reset(nullptr);
            }
            _executor.reset(nullptr);
            break;
        }
    }
}

const executor_base *const Config::executor() const {
    return _executor.get();
}

}
}
}