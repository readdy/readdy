/**
 * << detailed description >>
 *
 * @file SingleCPUEulerDBIntegrator.h
 * @brief << brief description >>
 * @author clonker
 * @date 19.04.16
 */

#ifndef READDY_MAIN_SingleCPUEulerBDIntegrator_H
#define READDY_MAIN_SingleCPUEulerBDIntegrator_H

#include <readdy/model/programs/Programs.h>
#include <readdy/kernel/singlecpu/SCPUKernel.h>

namespace readdy {
namespace kernel {
namespace scpu {

namespace programs {
class SCPUEulerBDIntegrator : public readdy::model::programs::EulerBDIntegrator {

public:
    SCPUEulerBDIntegrator(SCPUKernel *kernel);

    virtual void execute() override;

private:
    SCPUKernel *kernel;
};
}
}
}
}


#endif //READDY_MAIN_SingleCPUEulerBDIntegrator_H
