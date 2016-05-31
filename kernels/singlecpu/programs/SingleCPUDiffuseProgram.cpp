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
                    const auto& pos = pd->positions;
                    for (auto p = 0; p < pos->size(); p++) {
                        const double D = context.getDiffusionConstant((*pd->type)[p]);
                        auto displacement = sqrt(2. * D * dt) * (kernel->getRandomProvider().getNormal3());
                        (*pos)[p] += displacement;
                    }
                }
                SingleCPUDiffuseProgram::SingleCPUDiffuseProgram(SingleCPUKernel *kernel) : DiffuseProgram(), kernel(kernel) {};
            }
        }
    }
}


