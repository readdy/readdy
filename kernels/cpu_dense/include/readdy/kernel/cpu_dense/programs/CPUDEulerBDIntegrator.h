/**
 * << detailed description >>
 *
 * @file CPUEulerBDIntegrator.h
 * @brief << brief description >>
 * @author clonker
 * @date 22.11.16
 */

#ifndef READDY_DENSE_CPUEULERBDINTEGRATOR_H
#define READDY_DENSE_CPUEULERBDINTEGRATOR_H

#include <readdy/kernel/cpu_dense/CPUDKernel.h>
#include <readdy/model/programs/Programs.h>

namespace readdy {
namespace kernel {
namespace cpu_dense {
namespace programs {
class CPUDEulerBDIntegrator : public readdy::model::programs::EulerBDIntegrator {

public:
    CPUDEulerBDIntegrator(CPUDKernel *kernel);

    virtual void execute() override;

private:
    CPUDKernel *kernel;
};
}
}
}
}
#endif //READDY_DENSE_CPUEULERBDINTEGRATOR_H
