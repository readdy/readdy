/**
 * << detailed description >>
 *
 * @file SingleCPUDiffuseProgram.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 19.04.16
 */

#include "SingleCPUDiffuseProgram.h"
#include <boost/predef.h>

#if BOOST_OS_MACOS
#include <math.h>
#endif

void readdy::kernel::singlecpu::programs::SingleCPUDiffuseProgram::execute() {
    const auto& context = kernel->getKernelContext();
    const auto&& dt = context.getTimeStep();
    const auto&& pd = kernel->getKernelStateModelSingleCPU().getParticleData();
    const auto& pos = pd->positions;
    for (auto p = 0; p < pos->size(); p++) {
        const double D = context.getDiffusionConstant((*pd->type)[p]);
        const double prefactor = sqrt(2. * D * dt);

        auto displacement = kernel->getRandomProvider().getNormal3();
        displacement *= prefactor;
        (*pos)[p] += displacement;

    }
}

readdy::kernel::singlecpu::programs::SingleCPUDiffuseProgram::SingleCPUDiffuseProgram(SingleCPUKernel *kernel) : DiffuseProgram(), kernel(kernel) {};


