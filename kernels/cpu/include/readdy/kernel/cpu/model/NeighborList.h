/**
 * << detailed description >>
 *
 * @file NeighborList.h
 * @brief << brief description >>
 * @author clonker
 * @date 14.07.16
 */

#ifndef READDY_MAIN_NEIGHBORLIST_H
#define READDY_MAIN_NEIGHBORLIST_H

#include <memory>
#include <readdy/common/make_unique.h>
#include <readdy/kernel/cpu/model/ParticleIndexPair.h>
#include <readdy/model/KernelContext.h>
#include <readdy/kernel/singlecpu/model/SingleCPUNeighborList.h>

namespace readdy {
    namespace kernel {
        namespace cpu {
            namespace model {
                class NeighborList : public readdy::kernel::singlecpu::model::NotThatNaiveSingleCPUNeighborList<std::vector<ParticleIndexPair>> {

                public:
                    NeighborList(const readdy::model::KernelContext *const context) : NotThatNaiveSingleCPUNeighborList(context) { }
                };
            }
        }
    }
}
#endif //READDY_MAIN_NEIGHBORLIST_H
