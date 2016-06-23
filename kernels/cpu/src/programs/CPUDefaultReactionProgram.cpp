/**
 * << detailed description >>
 *
 * @file CPUDefaultReactionProgram.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 23.06.16
 */

#include <readdy/kernel/cpu/programs/CPUDefaultReactionProgram.h>

using super = readdy::kernel::singlecpu::programs::SingleCPUDefaultReactionProgram;

namespace readdy {
    namespace kernel {
        namespace cpu {
            namespace programs {
                CPUDefaultReactionProgram::CPUDefaultReactionProgram(const CPUKernel *const kernel) : super::SingleCPUDefaultReactionProgram(kernel)
                {

                }

                void CPUDefaultReactionProgram::execute() {
                    BOOST_LOG_TRIVIAL(debug) << "calling CPU default reaction program...";
                    super::execute();
                }

            }
        }
    }
}