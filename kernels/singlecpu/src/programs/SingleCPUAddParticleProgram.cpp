/**
 * << detailed description >>
 *
 * @file SingleCPUAddParticleProgram.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 19.04.16
 */

#include <readdy/kernel/singlecpu/programs/SingleCPUAddParticleProgram.h>

namespace kern = readdy::kernel::singlecpu::programs;

void kern::SingleCPUAddParticleProgram::execute() {
    kernel->getKernelStateModel().addParticles(particles);
}

readdy::kernel::singlecpu::programs::SingleCPUAddParticleProgram::SingleCPUAddParticleProgram(SingleCPUKernel *kernel)
        : kernel(kernel) {
}


kern::SingleCPUAddParticleProgram::~SingleCPUAddParticleProgram() = default;





