/**
 * << detailed description >>
 *
 * @file ConfigUtils.h
 * @brief << brief description >>
 * @author clonker
 * @date 01.08.16
 */

#ifndef READDY_MAIN_CONFIGUTILS_H
#define READDY_MAIN_CONFIGUTILS_H

#include <thread>

namespace readdy {
    namespace kernel {
        namespace cpu {
            namespace util {
                static unsigned long nThreads = std::thread::hardware_concurrency();

                inline unsigned long getNThreads() {
                    return nThreads;
                }

                inline void setNThreads(unsigned long threads) {
                    nThreads = threads < std::thread::hardware_concurrency() ?
                               threads : std::thread::hardware_concurrency();
                }
            }
        }
    }
}
#endif //READDY_MAIN_CONFIGUTILS_H
