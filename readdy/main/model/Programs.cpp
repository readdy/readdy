/**
 * << detailed description >>
 *
 * @file Programs.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 17.11.16
 */

#include <readdy/model/programs/Programs.h>
#include <readdy/common/logging.h>

namespace readdy {
namespace model {
namespace programs {


void UpdateNeighborList::setSkinSize(double skinSize) {
    log::console()->warn("The selected kernel has no Verlet list implemented, thus ignoring the skin size");
}
}
}
}