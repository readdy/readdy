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

#include "readdy/kernel/cpu/util/Config.h"

#if READDY_OSX
#include <cstdlib>
#endif

namespace readdy {
namespace kernel {
namespace cpu {
namespace util {
Config::Config() {
    nThreads = std::thread::hardware_concurrency();
    const char *env = std::getenv("READDY_N_CORES");
    if (env) {
        nThreads = static_cast<n_threads_t>(std::stol(env));
        log::console()->debug("Using {} threads (by environment variable READDY_N_CORES", nThreads);
    }
    nThreads = 1;
}
}
}
}
}