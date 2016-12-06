/**
 * << detailed description >>
 *
 * @file Gillespie.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 23.11.16
 */

#include <readdy/kernel/cpu_dense/programs/reactions/CPUDGillespie.h>


namespace readdy {
namespace kernel {
namespace cpu_dense {
namespace programs {
namespace reactions {

CPUDGillespie::CPUDGillespie(const CPUDKernel *const kernel) : kernel(kernel) {}

void CPUDGillespie::execute() {
    const auto &ctx = kernel->getKernelContext();
    auto data = kernel->getKernelStateModel().getParticleData();
    const auto &dist = ctx.getDistSquaredFun();
    const auto &fixPos = ctx.getFixPositionFun();
    const auto nl = kernel->getKernelStateModel().getNeighborList();

    double alpha = 0.0;
    std::vector<event_t> events;
    gatherEvents(kernel, readdy::util::range<event_t::index_type>(0, data->size()),
                 nl, *data, alpha, events, dist);
    auto particlesUpdate = handleEventsGillespie(kernel, false, true, std::move(events));

    // update data structure
    data->deactivateMarked();
    data->update(std::move(particlesUpdate));
}

}
}
}
}
}