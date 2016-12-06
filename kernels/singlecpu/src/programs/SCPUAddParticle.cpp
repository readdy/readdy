/**
 * << detailed description >>
 *
 * @file SingleCPUAddParticleProgram.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 19.04.16
 */

#include <readdy/kernel/singlecpu/programs/SCPUAddParticle.h>

namespace kern = readdy::kernel::scpu::programs;

void kern::SCPUAddParticle::execute() {
    kernel->getKernelStateModel().addParticles(particles);
}

readdy::kernel::scpu::programs::SCPUAddParticle::SCPUAddParticle(SCPUKernel *kernel)
        : kernel(kernel) {
}


kern::SCPUAddParticle::~SCPUAddParticle() = default;





