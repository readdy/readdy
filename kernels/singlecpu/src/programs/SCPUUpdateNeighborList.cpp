/**
 * << detailed description >>
 *
 * @file SingleCPUUpdateNeighborList.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 11.07.16
 */
#include <readdy/kernel/singlecpu/programs/SCPUUpdateNeighborList.h>

namespace readdy {
namespace kernel {
namespace scpu {
namespace programs {
SCPUUpdateNeighborList::SCPUUpdateNeighborList(SCPUKernel *kernel) : kernel(kernel) {}

void SCPUUpdateNeighborList::execute() {
    switch (action) {
        case create:
            kernel->getKernelStateModel().updateNeighborList();
            break;
        case clear:
            kernel->getKernelStateModel().clearNeighborList();
            break;
    }
}
}
}
}
}
