//
// Created by clonker on 07.03.16.
//

#ifndef READDY2_MAIN_SINGLECPUKERNEL_H
#define READDY2_MAIN_SINGLECPUKERNEL_H

#include <Kernel.h>

namespace readdy {
    namespace kernel {
        class SingleCPUKernel : readdy::plugin::Kernel {
            std::string getName();
        };
    }
}

#endif //READDY2_MAIN_SINGLECPUKERNEL_H
