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
    model->addParticles(particles);
}

readdy::kernel::singlecpu::programs::SingleCPUAddParticleProgram::SingleCPUAddParticleProgram(std::shared_ptr<SingleCPUKernelStateModel> model) {
    this->model = model;
}


k::SingleCPUAddParticleProgram::~SingleCPUAddParticleProgram() = default;





