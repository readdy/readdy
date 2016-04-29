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
#include <cmath>
#include <readdy/model/Vec3.h>

#if BOOST_OS_MACOS
#include <math.h>
#endif

void readdy::kernel::singlecpu::programs::SingleCPUDiffuseProgram::execute() {
    auto context = kernel->getKernelContext();
    auto randomProvider = kernel->getRandomProvider();
    auto dt = context->getTimeStep();
    auto pd = kernel->getKernelStateModelSingleCPU().getParticleData();
    auto pos = pd->positions;
    for (auto p = 0; p < pos->size(); p++) {
        const double D = context->getDiffusionConstant((*pd->type)[p]);
        const double prefactor = sqrt(2. * D * dt);

        readdy::model::Vec3 displacement {randomProvider->getNormal(), randomProvider->getNormal(), randomProvider->getNormal()};
        displacement *= prefactor;
        (*pos)[p] += displacement;

    }
}

readdy::kernel::singlecpu::programs::SingleCPUDiffuseProgram::SingleCPUDiffuseProgram(SingleCPUKernel *kernel) : DiffuseProgram(), kernel(kernel) {};


