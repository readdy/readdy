/**
 * << detailed description >>
 *
 * @file CPUEulerBDIntegrator.h
 * @brief << brief description >>
 * @author clonker
 * @date 07.07.16
 */

#ifndef READDY_MAIN_CPUEULERBDINTEGRATOR_H
#define READDY_MAIN_CPUEULERBDINTEGRATOR_H

#include <readdy/kernel/cpu/CPUKernel.h>
#include <readdy/model/programs/Programs.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace programs {
class CPUEulerBDIntegrator : public readdy::model::programs::EulerBDIntegrator {

public:
    CPUEulerBDIntegrator(CPUKernel *kernel);

    virtual void execute() override;

private:
    CPUKernel *kernel;
};
}
}
}
}
#endif //READDY_MAIN_CPUEULERBDINTEGRATOR_H
