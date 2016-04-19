//
// Created by clonker on 11.04.16.
//

#ifndef READDY2_MAIN_SINGLECPUTESTPROGRAM_H
#define READDY2_MAIN_SINGLECPUTESTPROGRAM_H

#include <readdy/plugin/Programs.h>

namespace readdy {
    namespace kernel {
        namespace singlecpu {
            namespace programs {
                class SingleCPUTestProgram : public readdy::plugin::TestProgram {
                public:
                    SingleCPUTestProgram();
                    ~SingleCPUTestProgram();

                    virtual void execute() override;
                };
            }
        }
    }
}

#endif //READDY2_MAIN_SINGLECPUTESTPROGRAM_H
