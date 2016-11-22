/**
 * << detailed description >>
 *
 * @file Gillespie.h
 * @brief << brief description >>
 * @author clonker
 * @date 20.10.16
 */

#ifndef READDY_DENSE_GILLESPIE_H
#define READDY_DENSE_GILLESPIE_H

#include <readdy/kernel/cpu_dense/Kernel.h>
#include <readdy/common/range.h>
#include "ReactionUtils.h"

namespace readdy {
namespace kernel {
namespace cpu_dense {
namespace programs {
namespace reactions {

class Gillespie : public readdy::model::programs::reactions::Gillespie {
    using event_t = Event;
    using reaction_idx_t = event_t::index_type;

public:

    Gillespie(Kernel const *const kernel);

    virtual void execute() override {
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

protected:
    Kernel const *const kernel;
};
}
}
}
}
}

#endif //READDY_DENSE_GILLESPIE_H
