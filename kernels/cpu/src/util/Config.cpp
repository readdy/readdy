/**
 * << detailed description >>
 *
 * @file Config.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 05.09.16
 */

#include <thread>
#include "readdy/kernel/cpu/util/Config.h"

namespace readdy {
    namespace kernel {
        namespace cpu {
            namespace util {
                Config::Config() {
                    nThreads = std::thread::hardware_concurrency();
                }
            }
        }
    }
}