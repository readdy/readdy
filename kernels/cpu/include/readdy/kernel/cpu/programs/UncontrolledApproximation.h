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
#include <readdy/kernel/singlecpu/programs/SingleCPUReactionImpls.h>

namespace readdy {
    namespace kernel {
        namespace cpu {
            namespace programs {
                namespace reactions {
                    class UncontrolledApproximation : public singlecpu::programs::reactions::UncontrolledApproximation {
                    public:
                        UncontrolledApproximation(const CPUKernel *const kernel);

                        virtual void execute() override;

                    protected:
                        CPUKernel const *const kernel;
                    };
                }

            }
        }
    }
}
#endif //READDY_MAIN_CPUDEFAULTREACTIONPROGRAM_H
