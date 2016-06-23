/**
 * << detailed description >>
 *
 * @file CPUDefaultReactionProgram.h
 * @brief << brief description >>
 * @author clonker
 * @date 23.06.16
 */

#ifndef READDY_MAIN_CPUDEFAULTREACTIONPROGRAM_H
#define READDY_MAIN_CPUDEFAULTREACTIONPROGRAM_H

#include <readdy/kernel/cpu/CPUKernel.h>
#include <readdy/kernel/singlecpu/programs/SingleCPUDefaultReactionProgram.h>

namespace readdy {
    namespace kernel {
        namespace cpu {
            namespace programs {
                class CPUDefaultReactionProgram : public singlecpu::programs::SingleCPUDefaultReactionProgram {
                public:
                    CPUDefaultReactionProgram(const CPUKernel *const kernel);

                    virtual void execute() override;

                };

            }
        }
    }
}
#endif //READDY_MAIN_CPUDEFAULTREACTIONPROGRAM_H
