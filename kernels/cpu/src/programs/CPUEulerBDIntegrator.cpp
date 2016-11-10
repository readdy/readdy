/**
 * << detailed description >>
 *
 * @file CPUEulerBDIntegrator.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 07.07.16
 */

#include <readdy/kernel/cpu/programs/CPUEulerBDIntegrator.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace programs {

namespace rnd = readdy::model::rnd;

void CPUEulerBDIntegrator::execute() {
    const auto &&pd = kernel->getKernelStateModel().getParticleData();
    const auto size = pd->size();
    std::vector<util::scoped_thread> threads;
    threads.reserve(kernel->getNThreads());
    const std::size_t grainSize = size / kernel->getNThreads();

    const auto &context = kernel->getKernelContext();
    using iter_t = decltype(pd->entries.begin());

    auto worker = [&context](iter_t entry_begin, iter_t entry_end)  {
        const auto &fixPos = context.getFixPositionFun();
        const auto kbt = context.getKBT();
        const auto dt = context.getTimeStep();
        for (auto it = entry_begin; it != entry_end; ++it) {
            if(!(*it).is_deactivated()) {
                const double D = context.getDiffusionConstant((*it).type);
                const auto randomDisplacement = sqrt(2. * D * dt) * rnd::normal3(0, 1);
                const auto deterministicDisplacement = (*it).force * dt * D / kbt;
                it->pos += randomDisplacement + deterministicDisplacement;
                fixPos(it->pos);
            }
        }

    };
    auto work_iter = pd->entries.begin();
    {
        for (unsigned int i = 0; i < kernel->getNThreads() - 1; ++i) {
            threads.push_back(util::scoped_thread(std::thread(worker, work_iter, work_iter + grainSize)));
            work_iter += grainSize;
        }
        threads.push_back(util::scoped_thread(std::thread(worker, work_iter, pd->entries.end())));
    }

}

CPUEulerBDIntegrator::CPUEulerBDIntegrator(CPUKernel *kernel) : kernel(kernel) {}

}
}
}
}