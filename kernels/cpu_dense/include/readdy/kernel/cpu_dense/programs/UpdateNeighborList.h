/**
 * << detailed description >>
 *
 * @file UpdateNeighborList.h
 * @brief << brief description >>
 * @author clonker
 * @date 22.11.16
 */


#ifndef READDY_DENSE_UPDATENEIGHBORLIST_H
#define READDY_DENSE_UPDATENEIGHBORLIST_H

#include <readdy/model/programs/Programs.h>
#include <readdy/kernel/cpu_dense/Kernel.h>

namespace readdy {
namespace kernel {
namespace cpu_dense {
namespace programs {
class UpdateNeighborList : public readdy::model::programs::UpdateNeighborList {
public:

    UpdateNeighborList(Kernel *kernel) : kernel(kernel) {
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
    Kernel *kernel;
};
}
}
}
}
#endif //READDY_DENSE_UPDATENEIGHBORLIST_H
