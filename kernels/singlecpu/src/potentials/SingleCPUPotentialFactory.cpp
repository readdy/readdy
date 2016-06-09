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
#include <readdy/kernel/singlecpu/potentials/PotentialsOrder2.h>

namespace readdy {
    namespace kernel {
        namespace singlecpu {
            namespace potentials {
                SingleCPUPotentialFactory::SingleCPUPotentialFactory(SingleCPUKernel *const kernel) : readdy::model::potentials::PotentialFactory(), kernel(kernel) {
                    factory[readdy::model::potentials::_internal::getPotentialName<HarmonicRepulsion>()] = [kernel] {return new HarmonicRepulsion(kernel);};
                }

            }
        }
    }
}