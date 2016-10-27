/**
 * << detailed description >>
 *
 * @file UncontrolledApproximation.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 20.10.16
 */

#include <readdy/kernel/cpu/programs/reactions/UncontrolledApproximation.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace programs {
namespace reactions {

using entry_type = readdy::kernel::cpu::model::ParticleData::Entry;

UncontrolledApproximation::UncontrolledApproximation(const CPUKernel *const kernel)
        : kernel(kernel) {

}

void UncontrolledApproximation::execute() {
    const auto &ctx = kernel->getKernelContext();
    const auto &fixPos = ctx.getFixPositionFun();
    const auto &dt = ctx.getTimeStep();
    auto data = kernel->getKernelStateModel().getParticleData();

    std::vector<entry_type> newParticles{};
    std::vector<std::function<void()>> events{};

    // shuffle reactions
    std::random_shuffle(events.begin(), events.end());

    // execute reactions
    std::for_each(events.begin(), events.end(), [](const std::function<void()> &f) { f(); });

    /**
     * TODO: This needs to be reworked so that it can function with inplace reactions. Also, efficiency.
     */

    // reposition particles to respect the periodic b.c.
    /*std::for_each(newParticles.begin(), newParticles.end(),
                  [&fixPos](particle_type &p) { fixPos(p.getPos()); });*/

    // update data structure
    /*data->addParticles(newParticles);*/
}
}
}
}
}
}