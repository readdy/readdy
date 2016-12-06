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
#include <readdy/kernel/singlecpu/SCPUKernel.h>

namespace readdy {
namespace kernel {
namespace scpu {
namespace programs {

class SCPUUpdateNeighborList : public readdy::model::programs::UpdateNeighborList {

public:
    SCPUUpdateNeighborList(SCPUKernel *kernel);

    virtual void execute() override;

private:
    SCPUKernel *kernel;
};

}
}
}
}

#endif //READDY_MAIN_SINGLECPUUPDATENEIGHBORLIST_H
