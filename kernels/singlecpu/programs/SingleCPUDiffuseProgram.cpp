/**
 * << detailed description >>
 *
 * @file SingleCPUDiffuseProgram.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 19.04.16
 */

#include "SingleCPUDiffuseProgram.h"
#include <boost/predef.h>

#if BOOST_OS_MACOS
#include <math.h>
#endif

void readdy::kernel::singlecpu::programs::SingleCPUDiffuseProgram::execute() {
    auto dt = context->getTimeStep();
    auto pd = model->particleData;
    auto pos = pd->positions;
    for (size_t p = 0; p < pos->size(); p++) {
        const double D = context->getDiffusionConstant((*pd->type)[p]);
        const double prefactor = sqrt(2. * D * dt);

        readdy::model::Vec3 displacement {randomProvider->getNormal(), randomProvider->getNormal(), randomProvider->getNormal()};
        displacement *= prefactor;
        (*pos)[p] += displacement;

    }
}

readdy::kernel::singlecpu::programs::SingleCPUDiffuseProgram::SingleCPUDiffuseProgram(std::shared_ptr<readdy::model::KernelContext> context, std::shared_ptr<SingleCPUKernelStateModel> model,
                                                                                      std::shared_ptr<readdy::utils::RandomProvider> randomProvider) : DiffuseProgram() {
    this->context = context;
    this->model = model;
    this->randomProvider = randomProvider;
}



