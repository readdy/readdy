/**
 * << detailed description >>
 *
 * @file Config.h
 * @brief << brief description >>
 * @author clonker
 * @date 05.09.16
 */

#ifndef READDY_CPUKERNEL_CONFIG_H
#define READDY_CPUKERNEL_CONFIG_H

namespace readdy {
namespace kernel {
namespace cpu {
namespace util {
struct Config {
    Config();

    unsigned long nThreads;
};
}
}
}
}
#endif //READDY_CPUKERNEL_CONFIG_H
