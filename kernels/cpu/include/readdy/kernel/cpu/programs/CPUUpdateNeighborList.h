/**
 * << detailed description >>
 *
 * @file UpdateNeighborList.h
 * @brief << brief description >>
 * @author clonker
 * @date 13.07.16
 */


#ifndef READDY_CPUKERNEL_UPDATENEIGHBORLIST_H
#define READDY_CPUKERNEL_UPDATENEIGHBORLIST_H

#include <readdy/model/programs/Programs.h>
#include <readdy/kernel/cpu/CPUKernel.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace programs {
class CPUUpdateNeighborList : public readdy::model::programs::UpdateNeighborList {
public:

    CPUUpdateNeighborList(CPUKernel *kernel) : kernel(kernel) {
    }

    virtual void execute() override {
        switch (action) {
            case create:
                kernel->getKernelStateModel().updateNeighborList();
                break;
            case clear:
                kernel->getKernelStateModel().clearNeighborList();
                break;
        }

    }

    virtual void setSkinSize(double skinSize) override {
        kernel->getKernelStateModel().getNeighborList()->setSkinSize(skinSize);
    }

private:
    CPUKernel *kernel;
};
}
}
}
}
#endif //READDY_CPUKERNEL_UPDATENEIGHBORLIST_H
