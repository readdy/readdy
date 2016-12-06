/**
 * << detailed description >>
 *
 * @file SingleCPUCalculateForces.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 11.07.16
 */

#include <readdy/kernel/singlecpu/programs/SCPUCalculateForces.h>

namespace readdy {
namespace kernel {
namespace scpu {
namespace programs {
SCPUCalculateForces::SCPUCalculateForces(SCPUKernel *kernel) : kernel(kernel) {}

void SCPUCalculateForces::execute() {
    kernel->getKernelStateModel().calculateForces();
}

}
}
}
}