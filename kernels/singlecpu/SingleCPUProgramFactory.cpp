//
// Created by clonker on 08.04.16.
//

#include "SingleCPUProgramFactory.h"
#include "programs/SingleCPUTestProgram.h"
#include "programs/SingleCPUAddParticleProgram.h"
#include "programs/SingleCPUDiffuseProgram.h"

std::shared_ptr<readdy::plugin::Program> readdy::kernel::singlecpu::SingleCPUProgramFactory::createProgram(const std::string name) {
    namespace prog = readdy::kernel::singlecpu::programs;
    if(name == prog::SingleCPUTestProgram::getName()) {
        return std::make_shared<prog::SingleCPUTestProgram>();
    }
    if(name == prog::SingleCPUAddParticleProgram::getName()) {
        return std::make_shared<prog::SingleCPUAddParticleProgram>(model);
    }
    if(name == prog::SingleCPUDiffuseProgram::getName()) {
        return std::make_shared<prog::SingleCPUDiffuseProgram>(context, model, randomProvider);
    }
    return nullptr;
}

readdy::kernel::singlecpu::SingleCPUProgramFactory::SingleCPUProgramFactory(std::shared_ptr<readdy::model::KernelContext> context, std::shared_ptr<SingleCPUKernelStateModel> model, std::shared_ptr<readdy::utils::RandomProvider> randomProvider) {
    this->context = context;
    this->model = model;
    this->randomProvider = randomProvider;
}


