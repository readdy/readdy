/**
 * << detailed description >>
 *
 * @file SingleCPUNeighborList.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 09.06.16
 */

#include <readdy/kernel/singlecpu/SCPUKernel.h>

namespace readdy {
namespace kernel {
namespace scpu {
namespace model {

void SCPUNaiveNeighborList::create(const SCPUParticleData &data) {
    pairs->clear();
    for (size_t i = 0; i < data.size(); ++i) {
        for (size_t j = i + 1; j < data.size(); ++j) {
            pairs->emplace(i, j);
        }
    }
}

SCPUNaiveNeighborList &SCPUNaiveNeighborList::operator=(
        SCPUNaiveNeighborList &&rhs) = default;


SCPUNaiveNeighborList::SCPUNaiveNeighborList() = default;

SCPUNaiveNeighborList::~SCPUNaiveNeighborList() = default;

SCPUNaiveNeighborList::SCPUNaiveNeighborList(SCPUNaiveNeighborList &&rhs) = default;
}
}
}
}
