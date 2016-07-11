/**
 * << detailed description >>
 *
 * @file CPUProgramFactory.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 23.06.16
 */

#include <readdy/kernel/cpu/programs/CPUProgramFactory.h>
#include <readdy/model/programs/Programs.h>
#include <readdy/kernel/cpu/programs/CPUDefaultReactionProgram.h>
#include <readdy/kernel/cpu/programs/CPUEulerBDIntegrator.h>

using super = readdy::kernel::singlecpu::programs::SingleCPUProgramFactory;
namespace core_p = readdy::model::programs;

namespace readdy {
    namespace kernel {
        namespace cpu {
            namespace programs {
                CPUProgramFactory::CPUProgramFactory(CPUKernel *kernel) : super::SingleCPUProgramFactory(kernel) {
                    factory[core_p::getProgramName<core_p::reactions::UncontrolledApproximation>()] = [kernel] {
                        return new CPUDefaultReactionProgram(kernel);
                    };
                    factory[core_p::getProgramName<core_p::EulerBDIntegrator>()] = [kernel] {
                        return new CPUEulerBDIntegrator(kernel);
                    };
                }
            }
        }
    }
}