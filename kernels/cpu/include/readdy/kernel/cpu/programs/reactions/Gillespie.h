/**
 * << detailed description >>
 *
 * @file Gillespie.h
 * @brief << brief description >>
 * @author clonker
 * @date 20.10.16
 */

#ifndef READDY_CPUKERNEL_GILLESPIE_H
#define READDY_CPUKERNEL_GILLESPIE_H

#include <readdy/kernel/cpu/CPUKernel.h>
#include <readdy/kernel/singlecpu/programs/SingleCPUReactionImpls.h>
#include <readdy/common/range.h>
#include "ReactionUtils.h"

namespace readdy {
namespace kernel {
namespace cpu {
namespace programs {
namespace reactions {

class Gillespie : public readdy::model::programs::reactions::Gillespie {
    using event_t = readdy::kernel::singlecpu::programs::reactions::ReactionEvent;
    using reaction_idx_t = event_t::index_type;

public:

    Gillespie(CPUKernel const *const kernel);

    virtual void execute() override {
        const auto &ctx = kernel->getKernelContext();
        auto data = kernel->getKernelStateModel().getParticleData();
        const auto &dist = ctx.getDistSquaredFun();
        const auto &fixPos = ctx.getFixPositionFun();
        const auto nl = kernel->getKernelStateModel().getNeighborList();

        double alpha = 0.0;
        std::vector<event_t> events;
        gatherEvents(kernel, readdy::util::range<std::size_t>(0, data->size()), nl, *data, alpha, events);
        auto newParticles = handleEventsGillespie(kernel, false, true, std::move(events));

        // reposition particles to respect the periodic b.c.
        std::for_each(newParticles.begin(), newParticles.end(),
                      [&fixPos](data_t::Entry &p) { fixPos(p.pos); });

        // update data structure
        data->addEntries(newParticles);
    }

protected:
    CPUKernel const *const kernel;
};
}
}
}
}
}

#endif //READDY_CPUKERNEL_GILLESPIE_H
