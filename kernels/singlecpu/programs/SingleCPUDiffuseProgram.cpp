/**
 * << detailed description >>
 *
 * @file SingleCPUDiffuseProgram.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 19.04.16
 */

#include "SingleCPUDiffuseProgram.h"
#include "../SingleCPUKernelStateModel.h"

readdy::kernel::singlecpu::programs::SingleCPUDiffuseProgram::SingleCPUDiffuseProgram(readdy::kernel::singlecpu::SingleCPUKernel &kernel) : DiffuseProgram(){
    SingleCPUDiffuseProgram::kernel = kernel;
}

void readdy::kernel::singlecpu::programs::SingleCPUDiffuseProgram::execute() {
    auto rand = kernel.getRandomProvider();
    auto kernelModel = dynamic_cast<readdy::kernel::singlecpu::SingleCPUKernelStateModel*>(kernel.getKernelStateModel().get());
    auto ctx = kernel.getKernelContext();
    auto dt = ctx->getTimeStep();
    for (auto &&particle :kernelModel->getParticles()) {
        const double D = ctx->getDiffusionConstant(particle.type);
        const double prefactor = sqrt(2. * D * dt);
        particle.pos[0] += prefactor * rand->getNormal();
        particle.pos[1] += prefactor * rand->getNormal();
        particle.pos[2] += prefactor * rand->getNormal();
    }
}

