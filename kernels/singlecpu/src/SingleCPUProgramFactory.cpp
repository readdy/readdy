//
// Created by clonker on 08.04.16.
//

#include <readdy/common/make_unique.h>
#include <readdy/kernel/singlecpu/SingleCPUProgramFactory.h>
#include <readdy/kernel/singlecpu/programs/SingleCPUTestProgram.h>
#include <readdy/kernel/singlecpu/programs/SingleCPUAddParticleProgram.h>
#include <readdy/kernel/singlecpu/programs/SingleCPUDiffuseProgram.h>
#include <readdy/kernel/singlecpu/programs/SingleCPUUpdateStateModelProgram.h>

std::unique_ptr<readdy::model::Program> readdy::kernel::singlecpu::SingleCPUProgramFactory::createProgram(const std::string& name) const {
    namespace prog = readdy::kernel::singlecpu::programs;
    namespace core_p = readdy::model;
    if(name == core_p::getProgramName<core_p::TestProgram>()) {
        return std::make_unique<prog::SingleCPUTestProgram>();
    }
    if(name == core_p::getProgramName<core_p::AddParticleProgram>()) {
        return std::make_unique<prog::SingleCPUAddParticleProgram>(kernel);
    }
    if(name == core_p::getProgramName<core_p::DiffuseProgram>()) {
        return std::make_unique<prog::SingleCPUDiffuseProgram>(kernel);
    }
    if(name == core_p::getProgramName<core_p::UpdateStateModelProgram>()) {
        return std::make_unique<prog::SingleCPUUpdateStateModelProgram>(kernel);
    }
    return nullptr;
}

readdy::kernel::singlecpu::SingleCPUProgramFactory::SingleCPUProgramFactory(SingleCPUKernel *kernel) : kernel(kernel){

}


