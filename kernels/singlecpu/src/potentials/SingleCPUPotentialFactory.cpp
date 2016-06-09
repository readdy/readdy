/**
 * << detailed description >>
 *
 * @file SingleCPUPotentialFactory.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 09.06.16
 */

#include <readdy/kernel/singlecpu/SingleCPUKernel.h>
#include <readdy/kernel/singlecpu/potentials/SingleCPUPotentialFactory.h>

namespace readdy {
    namespace kernel {
        namespace singlecpu {
            namespace potentials {
                SingleCPUPotentialFactory::SingleCPUPotentialFactory(SingleCPUKernel *const kernel) : readdy::model::potentials::PotentialFactory(), kernel(kernel) {
                    // todo add potentials
                }

            }
        }
    }
}