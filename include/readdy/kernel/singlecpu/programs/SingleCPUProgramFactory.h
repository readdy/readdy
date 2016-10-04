//
// Created by clonker on 08.04.16.
//

#ifndef READDY_MAIN_SINGLECPUPROGRAMFACTORY_H
#define READDY_MAIN_SINGLECPUPROGRAMFACTORY_H

#include <readdy/model/programs/ProgramFactory.h>
#include <readdy/kernel/singlecpu/SingleCPUKernelStateModel.h>

namespace readdy {
namespace kernel {
namespace singlecpu {
class SingleCPUKernel;
namespace programs {
class SingleCPUProgramFactory : public readdy::model::programs::ProgramFactory {
public:
    SingleCPUProgramFactory(SingleCPUKernel *kernel);

private:
    SingleCPUKernel *kernel;
};
}
}
}
}


#endif //READDY_MAIN_SINGLECPUPROGRAMFACTORY_H
