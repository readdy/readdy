/**
 * << detailed description >>
 *
 * @file SingleCPUEulerDBIntegrator.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 19.04.16
 */

#include <readdy/kernel/singlecpu/programs/SingleCPUEulerBDIntegrator.h>

#if BOOST_OS_MACOS
#include <math.h>
#endif
namespace readdy {
namespace kernel {
namespace singlecpu {
namespace programs {

SingleCPUEulerBDIntegrator::SingleCPUEulerBDIntegrator(SingleCPUKernel *kernel)
        : readdy::model::programs::EulerBDIntegrator(), kernel(kernel) {};

void SingleCPUEulerBDIntegrator::execute() {
    const auto &context = kernel->getKernelContext();
    const auto &kbt = context.getKBT();
    const auto &fixPos = context.getFixPositionFun();
    const auto &&dt = context.getTimeStep();
    const auto &&pd = kernel->getKernelStateModel().getParticleData();
    auto it_pos = pd->begin_positions();
    auto it_types = pd->begin_types();
    auto it_forces = pd->begin_forces();
    readdy::model::RandomProvider rnd;
    for (; it_pos != pd->end_positions();) {
        const double D = context.getDiffusionConstant(*it_types);
        const auto randomDisplacement = sqrt(2. * D * dt) * (rnd.getNormal3());
        *it_pos += randomDisplacement;
        const auto deterministicDisplacement = *it_forces * dt * D / kbt;
        *it_pos += deterministicDisplacement;
        fixPos(*it_pos);
        ++it_pos;
        ++it_types;
        ++it_forces;
    }
}
}
}
}
}


