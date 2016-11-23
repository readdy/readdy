/**
 * << detailed description >>
 *
 * @file Gillespie.h
 * @brief << brief description >>
 * @author clonker
 * @date 22.11.16
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

    virtual void execute() override;

protected:
    Kernel const *const kernel;
};
}
}
}
}
}

#endif //READDY_DENSE_GILLESPIE_H
