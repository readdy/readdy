/**
 * << detailed description >>
 *
 * @file SingleCPUNeighborList.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 09.06.16
 */

#include <readdy/kernel/singlecpu/SingleCPUKernel.h>

namespace readdy {
    namespace kernel {
        namespace singlecpu {
            namespace model {

                void NaiveSingleCPUNeighborList::create(const SingleCPUParticleData &data) {
                    pairs->clear();
                    for (size_t i = 0; i < data.size(); ++i) {
                        for (size_t j = i + 1; j < data.size(); ++j) {
                            pairs->emplace(i, j);
                        }
                    }
                }

                NaiveSingleCPUNeighborList &NaiveSingleCPUNeighborList::operator=(
                        NaiveSingleCPUNeighborList &&rhs) = default;


                NaiveSingleCPUNeighborList::NaiveSingleCPUNeighborList() = default;

                NaiveSingleCPUNeighborList::~NaiveSingleCPUNeighborList() = default;

                NaiveSingleCPUNeighborList::NaiveSingleCPUNeighborList(NaiveSingleCPUNeighborList &&rhs) = default;
            }
        }
    }
}
