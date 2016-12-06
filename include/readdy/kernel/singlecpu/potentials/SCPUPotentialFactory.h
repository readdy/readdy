/**
 * @file SingleCPUPotentialFactory.h
 * @brief This factory creates potentials for the singlecpu kernel.
 * @author clonker
 * @date 09.06.16
 */

#include <readdy/model/potentials/PotentialFactory.h>

#ifndef READDY_MAIN_SINGLECPUPOTENTIALFACTORY_H
#define READDY_MAIN_SINGLECPUPOTENTIALFACTORY_H

class SingleCPUKernel;

namespace readdy {
namespace kernel {
namespace scpu {
namespace potentials {
class SCPUPotentialFactory : public readdy::model::potentials::PotentialFactory {

public:
    SCPUPotentialFactory(SCPUKernel *const kernel);

protected:
    SCPUKernel const *kernel;
};
}
}
}
}
#endif //READDY_MAIN_SINGLECPUPOTENTIALFACTORY_H
