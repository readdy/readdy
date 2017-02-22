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
#include <readdy/common/logging.h>
#include <readdy/common/macros.h>

#include "readdy/common/thread/Config.h"

#if READDY_OSX
#include <cstdlib>
#endif

namespace readdy {
namespace util {
namespace thread {

Config::Config() {
    // magic number 4 to enable some load balancing
#ifdef READDY_DEBUG
    m_nThreads = std::thread::hardware_concurrency();
#else
    m_nThreads = 4*std::thread::hardware_concurrency();
#endif

    const char *env = std::getenv("READDY_N_CORES");
    if (env) {
        m_nThreads = static_cast<n_threads_t>(std::stol(env));
        log::debug("Using {} threads (by environment variable READDY_N_CORES", m_nThreads);
    }
}

Config::n_threads_t Config::nThreads() const {
    return m_nThreads;
}

void Config::setNThreads(const Config::n_threads_t n) {
    m_nThreads = n;
}

}
}
}