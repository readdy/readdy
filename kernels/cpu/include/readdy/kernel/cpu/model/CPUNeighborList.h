/**
 * << detailed description >>
 *
 * @file CPUNeighborList.h
 * @brief << brief description >>
 * @author clonker
 * @date 11.07.16
 */


#ifndef READDY_MAIN_CPUNEIGHBORLIST_H
#define READDY_MAIN_CPUNEIGHBORLIST_H

#include <readdy/kernel/singlecpu/model/SingleCPUNeighborList.h>

namespace readdy {
    namespace kernel {
        namespace cpu {
            namespace model {
                struct NotThatNaiveCPUNeighborList : public readdy::kernel::singlecpu::model::NotThatNaiveSingleCPUNeighborList {

                    NotThatNaiveCPUNeighborList(const readdy::model::KernelContext *const ctx);

                    virtual void fillBoxes(const readdy::kernel::singlecpu::model::SingleCPUParticleData &data) override;

                };
            }
        }
    }
}
#endif //READDY_MAIN_CPUNEIGHBORLIST_H
