/**
 * << detailed description >>
 *
 * @file SincleCPUUpdateStateModelProgram.h
 * @brief << brief description >>
 * @author clonker
 * @date 20.06.16
 */

#include <readdy/model/programs/Programs.h>
#include <readdy/kernel/singlecpu/SingleCPUKernel.h>

#ifndef READDY_MAIN_SINGLECPUUPDATESTATEMODELPROGRAM_H
#define READDY_MAIN_SINGLECPUUPDATESTATEMODELPROGRAM_H

namespace readdy {
    namespace kernel {
        namespace singlecpu {
            namespace programs {
                class SingleCPUUpdateStateModelProgram : public readdy::model::programs::UpdateStateModelProgram {
                public:
                    SingleCPUUpdateStateModelProgram(SingleCPUKernel *kernel);
                    virtual void execute() override;

                private:
                    SingleCPUKernel *kernel;
                    readdy::model::time_step_type t;
                };
            }
        }
    }
}

#endif //READDY_MAIN_SINGLECPUUPDATESTATEMODELPROGRAM_H
