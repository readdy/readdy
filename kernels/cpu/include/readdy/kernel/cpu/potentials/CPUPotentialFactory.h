/**
 * @file CPUPotentialFactory.h
 * @brief This factory creates potentials for the cpu kernel.
 * @author clonker
 * @date 13.07.16
 */


#ifndef READDY_CPUKERNEL_CPUPOTENTIALFACTORY_H
#define READDY_CPUKERNEL_CPUPOTENTIALFACTORY_H

#include <readdy/model/potentials/PotentialFactory.h>

namespace readdy {
namespace kernel {
namespace cpu {
class CPUKernel;
namespace potentials {
class CPUPotentialFactory : public readdy::model::potentials::PotentialFactory {
public:
    CPUPotentialFactory(CPUKernel *const kernel);
};
}
}
}
}
#endif //READDY_CPUKERNEL_CPUPOTENTIALFACTORY_H
