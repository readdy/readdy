/**
 * << detailed description >>
 *
 * @file SingleCPUAddParticleProgram.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 19.04.16
 */

#include "SingleCPUAddParticleProgram.h"

namespace k = readdy::kernel::singlecpu::programs;

void k::SingleCPUAddParticleProgram::execute() {
    kernel->getKernelStateModel().addParticles(particles);
}

readdy::kernel::singlecpu::programs::SingleCPUAddParticleProgram::SingleCPUAddParticleProgram(SingleCPUKernel* kernel) : kernel(kernel){
}


k::SingleCPUAddParticleProgram::~SingleCPUAddParticleProgram() = default;





