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

}

readdy::kernel::singlecpu::programs::SingleCPUAddParticleProgram::SingleCPUAddParticleProgram(const SingleCPUKernel &kernel) : kernel(kernel){

}


k::SingleCPUAddParticleProgram::~SingleCPUAddParticleProgram() = default;





