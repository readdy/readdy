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
#include <readdy/kernel/singlecpu/SingleCPUKernel.h>

namespace readdy {
namespace kernel {
namespace singlecpu {
namespace programs {
class SingleCPUCalculateForces : public readdy::model::programs::CalculateForces {
public:
    SingleCPUCalculateForces(SingleCPUKernel *kernel);

    virtual void execute() override;

private:
    SingleCPUKernel *kernel;
};
}
}
}
}

#endif //READDY_MAIN_SINGLECPUCALCULATEFORCES_H
