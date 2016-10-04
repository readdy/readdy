//
// Created by clonker on 11.04.16.
//

#ifndef READDY_MAIN_SINGLECPUTESTPROGRAM_H
#define READDY_MAIN_SINGLECPUTESTPROGRAM_H

#include <readdy/model/programs/Programs.h>

namespace readdy {
namespace kernel {
namespace singlecpu {
namespace programs {
class SingleCPUTestProgram : public readdy::model::programs::Test {
public:
    SingleCPUTestProgram();

    ~SingleCPUTestProgram();

    virtual void execute() override;
};
}
}
}
}

#endif //READDY_MAIN_SINGLECPUTESTPROGRAM_H
