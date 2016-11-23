/**
 * << detailed description >>
 *
 * @file CPUEulerBDIntegrator.h
 * @brief << brief description >>
 * @author clonker
 * @date 07.07.16
 */

#ifndef READDY_CPUKERNEL_CPUEULERBDINTEGRATOR_H
#define READDY_CPUKERNEL_CPUEULERBDINTEGRATOR_H

#include <readdy/kernel/cpu/Kernel.h>
#include <readdy/model/programs/Programs.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace programs {
class EulerBDIntegrator : public readdy::model::programs::EulerBDIntegrator {

public:
    EulerBDIntegrator(Kernel *kernel);

    virtual void execute() override;

private:
    Kernel *kernel;
};
}
}
}
}
#endif //READDY_CPUKERNEL_CPUEULERBDINTEGRATOR_H
