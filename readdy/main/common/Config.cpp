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
    m_nThreads = std::thread::hardware_concurrency();
    const char *env = std::getenv("READDY_N_CORES");
    if (env) {
        m_nThreads = static_cast<n_threads_t>(std::stol(env));
        log::console()->debug("Using {} threads (by environment variable READDY_N_CORES", m_nThreads);
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