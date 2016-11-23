/**
 * << detailed description >>
 *
 * @file CPUEulerBDIntegrator.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 07.07.16
 */

#include <readdy/common/thread/scoped_thread.h>
#include <readdy/kernel/cpu/programs/EulerBDIntegrator.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace programs {

namespace rnd = readdy::model::rnd;
namespace thd = readdy::util::thread;

void EulerBDIntegrator::execute() {
    auto& pd = *kernel->getKernelStateModel().getParticleData();
    const auto size = pd.size();
    std::vector<thd::scoped_thread> threads;
    threads.reserve(kernel->getNThreads());
    const std::size_t grainSize = size / kernel->getNThreads();

    const auto &context = kernel->getKernelContext();
    using iter_t = decltype(pd.begin());

    auto worker = [&context, &pd](iter_t entry_begin, iter_t entry_end)  {
        const auto &fixPos = context.getFixPositionFun();
        const auto kbt = context.getKBT();
        const auto dt = context.getTimeStep();
        for (iter_t it = entry_begin; it != entry_end; ++it) {
            if(!it->is_deactivated()) {
                const double D = context.getDiffusionConstant(it->type);
                const auto randomDisplacement = std::sqrt(2. * D * dt) * rnd::normal3(0, 1);
                const auto deterministicDisplacement = it->force * dt * D / kbt;
                pd.displace(*it, randomDisplacement + deterministicDisplacement);
            }
        }

    };
    auto work_iter = pd.begin();
    {
        for (unsigned int i = 0; i < kernel->getNThreads() - 1; ++i) {
            threads.push_back(thd::scoped_thread(std::thread(worker, work_iter, work_iter + grainSize)));
            work_iter += grainSize;
        }
        threads.push_back(thd::scoped_thread(std::thread(worker, work_iter, pd.end())));
    }

}

EulerBDIntegrator::EulerBDIntegrator(Kernel *kernel) : kernel(kernel) {}

}
}
}
}