/**
 * << detailed description >>
 *
 * @file Config.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 05.09.16
 */

#include <thread>
#include <boost/log/trivial.hpp>
#include "readdy/kernel/cpu/util/Config.h"

namespace readdy {
    namespace kernel {
        namespace cpu {
            namespace util {
                Config::Config() {
                    nThreads = std::thread::hardware_concurrency();
                    const char *env = std::getenv("READDY_N_CORES");
                    if(env) {
                        nThreads = static_cast<unsigned long>(std::stol(env));
                        BOOST_LOG_TRIVIAL(debug) << "Using "
                                                 << nThreads << " threads (by environment variable READDY_N_CORES)";
                    }
                }
            }
        }
    }
}