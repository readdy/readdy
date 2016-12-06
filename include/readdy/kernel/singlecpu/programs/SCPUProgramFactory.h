//
// Created by clonker on 08.04.16.
//

#ifndef READDY_MAIN_SINGLECPUPROGRAMFACTORY_H
#define READDY_MAIN_SINGLECPUPROGRAMFACTORY_H

#include <readdy/model/programs/ProgramFactory.h>
#include <readdy/kernel/singlecpu/SCPUStateModel.h>

namespace readdy {
namespace kernel {
namespace scpu {
class SCPUKernel;
namespace programs {
class SCPUProgramFactory : public readdy::model::programs::ProgramFactory {
public:
    SCPUProgramFactory(SCPUKernel *kernel);

private:
    SCPUKernel *kernel;
};
}
}
}
}


#endif //READDY_MAIN_SINGLECPUPROGRAMFACTORY_H
