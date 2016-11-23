/**
 * << detailed description >>
 *
 * @file Gillespie.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 20.10.16
 */
#include <readdy/kernel/cpu/programs/reactions/Gillespie.h>


namespace readdy {
namespace kernel {
namespace cpu {
namespace programs {
namespace reactions {

Gillespie::Gillespie(const Kernel *const kernel) : kernel(kernel) {}

void Gillespie::execute() {
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
    nl->updateData(std::move(particlesUpdate));
}

}
}
}
}
}