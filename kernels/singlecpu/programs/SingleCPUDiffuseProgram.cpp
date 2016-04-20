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

void readdy::kernel::singlecpu::programs::SingleCPUDiffuseProgram::execute() {
    auto dt = context->getTimeStep();
    for (auto &&particle : *model->getParticles()) {
        const double D = context->getDiffusionConstant(particle.type);
        const double prefactor = sqrt(2. * D * dt);
        particle.pos[0] += prefactor * randomProvider->getNormal();
        particle.pos[1] += prefactor * randomProvider->getNormal();
        particle.pos[2] += prefactor * randomProvider->getNormal();
    }
}

readdy::kernel::singlecpu::programs::SingleCPUDiffuseProgram::SingleCPUDiffuseProgram(std::shared_ptr<readdy::model::KernelContext> context, std::shared_ptr<SingleCPUKernelStateModel> model,
                                                                                      std::shared_ptr<readdy::utils::RandomProvider> randomProvider) : DiffuseProgram() {
    this->context = context;
    this->model = model;
    this->randomProvider = randomProvider;
}



