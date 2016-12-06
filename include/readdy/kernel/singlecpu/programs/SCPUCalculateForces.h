/**
 * << detailed description >>
 *
 * @file SincleCPUUpdateStateModelProgram.h
 * @brief << brief description >>
 * @author clonker
 * @date 20.06.16
 */

#ifndef READDY_MAIN_SINGLECPUCALCULATEFORCES_H
#define READDY_MAIN_SINGLECPUCALCULATEFORCES_H

#include <readdy/model/programs/Programs.h>
#include <readdy/kernel/singlecpu/SCPUKernel.h>

namespace readdy {
namespace kernel {
namespace scpu {
namespace programs {
class SCPUCalculateForces : public readdy::model::programs::CalculateForces {
public:
    SCPUCalculateForces(SCPUKernel *kernel);

    virtual void execute() override;

private:
    SCPUKernel *kernel;
};
}
}
}
}

#endif //READDY_MAIN_SINGLECPUCALCULATEFORCES_H
