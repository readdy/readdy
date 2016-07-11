/**
 * << detailed description >>
 *
 * @file SingleCPUCalculateForces.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 11.07.16
 */

#include <readdy/kernel/singlecpu/programs/SingleCPUCalculateForces.h>

namespace readdy {
    namespace kernel {
        namespace singlecpu {
            namespace programs {
                SingleCPUCalculateForces::SingleCPUCalculateForces(SingleCPUKernel *kernel) : kernel(kernel) {}

                void SingleCPUCalculateForces::execute() {
                    kernel->getKernelStateModelSingleCPU().calculateForces();
                }

            }
        }
    }
}