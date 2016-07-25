/**
 * << detailed description >>
 *
 * @file SingleCPUPotentialFactory.h
 * @brief << brief description >>
 * @author clonker
 * @date 09.06.16
 */

#include <readdy/model/potentials/PotentialFactory.h>

#ifndef READDY_MAIN_SINGLECPUPOTENTIALFACTORY_H
#define READDY_MAIN_SINGLECPUPOTENTIALFACTORY_H

class SingleCPUKernel;

namespace readdy {
    namespace kernel {
        namespace singlecpu {
            namespace potentials {
                class SingleCPUPotentialFactory : public readdy::model::potentials::PotentialFactory {

                public:
                    SingleCPUPotentialFactory(SingleCPUKernel *const kernel);

                protected:
                    SingleCPUKernel const* kernel;
                };
            }
        }
    }
}
#endif //READDY_MAIN_SINGLECPUPOTENTIALFACTORY_H
