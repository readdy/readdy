/**
 * << detailed description >>
 *
 * @file CalculateForces.h
 * @brief << brief description >>
 * @author clonker
 * @date 14.07.16
 */


#ifndef READDY_CPUKERNEL_CALCULATEFORCES_H
#define READDY_CPUKERNEL_CALCULATEFORCES_H

#include <readdy/model/programs/Programs.h>
#include <readdy/kernel/cpu/Kernel.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace programs {
class CalculateForces : public readdy::model::programs::CalculateForces {

public:

    CalculateForces(Kernel *kernel) : kernel(kernel) {}

    virtual void execute() override {
        kernel->getKernelStateModel().calculateForces();
    }

protected:
    Kernel *kernel;
};
}
}
}
}
#endif //READDY_CPUKERNEL_CALCULATEFORCES_H
