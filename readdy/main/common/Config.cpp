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

Config::Config() : Config(ThreadMode::std_thread) {}

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

Config::Config(ThreadMode mode) : _mode(mode), m_nThreads(readdy_default_n_threads()) {
    update();
}

}
}
}