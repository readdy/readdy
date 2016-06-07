//
// Created by clonker on 08.04.16.
//

#include <readdy/common/make_unique.h>
#include <readdy/kernel/singlecpu/SingleCPUProgramFactory.h>
#include <readdy/kernel/singlecpu/programs/SingleCPUTestProgram.h>
#include <readdy/kernel/singlecpu/programs/SingleCPUAddParticleProgram.h>
#include <readdy/kernel/singlecpu/programs/SingleCPUDiffuseProgram.h>

std::unique_ptr<readdy::model::Program> readdy::kernel::singlecpu::SingleCPUProgramFactory::createProgram(const std::string& name) const {
    namespace prog = readdy::kernel::singlecpu::programs;
    if(name == prog::SingleCPUTestProgram::getName()) {
        return std::make_unique<prog::SingleCPUTestProgram>();
    }
    if(name == prog::SingleCPUAddParticleProgram::getName()) {
        return std::make_unique<prog::SingleCPUAddParticleProgram>(kernel);
    }
    if(name == prog::SingleCPUDiffuseProgram::getName()) {
        return std::make_unique<prog::SingleCPUDiffuseProgram>(kernel);
    }
    return nullptr;
}

readdy::kernel::singlecpu::SingleCPUProgramFactory::SingleCPUProgramFactory(SingleCPUKernel *kernel) : kernel(kernel){

}


