/**
 * << detailed description >>
 *
 * @file SingleCPUUpdateNeighborList.h
 * @brief << brief description >>
 * @author clonker
 * @date 11.07.16
 */

#ifndef READDY_MAIN_SINGLECPUUPDATENEIGHBORLIST_H
#define READDY_MAIN_SINGLECPUUPDATENEIGHBORLIST_H

#include <readdy/model/programs/Programs.h>
#include <readdy/kernel/singlecpu/SingleCPUKernel.h>

namespace readdy {
    namespace kernel {
        namespace singlecpu {
            namespace programs {

                class SingleCPUUpdateNeighborList : public readdy::model::programs::UpdateNeighborList {
                    
                public:
                    SingleCPUUpdateNeighborList(SingleCPUKernel* kernel);
                    virtual void execute() override;

                private:
                    SingleCPUKernel * kernel;
                };

            }
        }
    }
}

#endif //READDY_MAIN_SINGLECPUUPDATENEIGHBORLIST_H
