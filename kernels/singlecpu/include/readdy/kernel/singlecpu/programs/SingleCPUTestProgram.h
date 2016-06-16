//
// Created by clonker on 11.04.16.
//

#ifndef READDY_MAIN_SINGLECPUTESTPROGRAM_H
#define READDY_MAIN_SINGLECPUTESTPROGRAM_H

#include <readdy/model/Programs.h>

namespace readdy {
    namespace kernel {
        namespace singlecpu {
            namespace programs {
                class SingleCPUTestProgram : public readdy::model::TestProgram {
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
