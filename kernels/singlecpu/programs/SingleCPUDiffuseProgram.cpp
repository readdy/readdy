/**
 * << detailed description >>
 *
 * @file SingleCPUDiffuseProgram.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 19.04.16
 */

#include "SingleCPUDiffuseProgram.h"

#if BOOST_OS_MACOS
#include <math.h>
#endif
namespace readdy {
    namespace kernel {
        namespace singlecpu {
            namespace programs {
                void SingleCPUDiffuseProgram::execute() {
                    const auto& context = kernel->getKernelContext();
                    const auto&& dt = context.getTimeStep();
                    const auto&& pd = kernel->getKernelStateModelSingleCPU().getParticleData();
                    auto it_pos = pd->begin_positions();
                    auto it_types = pd->begin_types();
                    for (; it_pos != pd->end_positions() ; ) {
                        const double D = context.getDiffusionConstant(*it_types);
                        auto displacement = sqrt(2. * D * dt) * (kernel->getRandomProvider().getNormal3());
                        *it_pos += displacement;
                        it_pos = std::next(it_pos);
                        it_types = std::next(it_types);
                    }
                }
                SingleCPUDiffuseProgram::SingleCPUDiffuseProgram(SingleCPUKernel *kernel) : DiffuseProgram(), kernel(kernel) {};
            }
        }
    }
}


