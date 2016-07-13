/**
 * << detailed description >>
 *
 * @file SingleCPUUpdateNeighborList.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 11.07.16
 */
#include <readdy/kernel/singlecpu/programs/SingleCPUUpdateNeighborList.h>

namespace readdy {
    namespace kernel {
        namespace singlecpu {
            namespace programs {
                SingleCPUUpdateNeighborList::SingleCPUUpdateNeighborList(SingleCPUKernel *kernel) : kernel(kernel) { }

                void SingleCPUUpdateNeighborList::execute() {
                    kernel->getKernelStateModel().updateNeighborList();
                }
            }
        }
    }
}
