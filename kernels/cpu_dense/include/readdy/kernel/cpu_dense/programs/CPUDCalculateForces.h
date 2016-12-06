/**
 * << detailed description >>
 *
 * @file CalculateForces.h
 * @brief << brief description >>
 * @author clonker
 * @date 22.11.16
 */


#ifndef READDY_DENSE_CALCULATEFORCES_H
#define READDY_DENSE_CALCULATEFORCES_H

#include <readdy/model/programs/Programs.h>
#include <readdy/kernel/cpu_dense/CPUDKernel.h>

namespace readdy {
namespace kernel {
namespace cpu_dense {
namespace programs {
class CPUDCalculateForces : public readdy::model::programs::CalculateForces {

public:

    CPUDCalculateForces(CPUDKernel *kernel) : kernel(kernel) {}

    virtual void execute() override {
        kernel->getKernelStateModel().calculateForces();
    }

protected:
    CPUDKernel *kernel;
};
}
}
}
}
#endif //READDY_DENSE_CALCULATEFORCES_H
