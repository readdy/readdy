/**
 * << detailed description >>
 *
 * @file UpdateNeighborList.h
 * @brief << brief description >>
 * @author clonker
 * @date 13.07.16
 */


#ifndef READDY_MAIN_UPDATENEIGHBORLIST_H
#define READDY_MAIN_UPDATENEIGHBORLIST_H

#include <readdy/model/programs/Programs.h>
#include <readdy/kernel/cpu/CPUKernel.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace programs {
class UpdateNeighborList : public readdy::model::programs::UpdateNeighborList {
public:

    UpdateNeighborList(CPUKernel *kernel) : kernel(kernel) {
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

private:
    CPUKernel *kernel;
};
}
}
}
}
#endif //READDY_MAIN_UPDATENEIGHBORLIST_H
