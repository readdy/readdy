/**
 * << detailed description >>
 *
 * @file FilteredGillespieParallel.h
 * @brief << brief description >>
 * @author clonker
 * @date 22.11.16
 */

#ifndef READDY_DENSE_FILTEREDGILLESPIEPARALLEL_H
#define READDY_DENSE_FILTEREDGILLESPIEPARALLEL_H

#include "GillespieParallel.h"

namespace readdy {
namespace kernel {
namespace cpu_dense {
namespace programs {
namespace reactions {

class FilteredGillespieParallel : public GillespieParallel {
public:
    FilteredGillespieParallel(const kernel_t *const kernel);

private:
    virtual void handleBoxReactions() override;
};

}
}
}
}
}
#endif //READDY_DENSE_FILTEREDGILLESPIEPARALLEL_H
