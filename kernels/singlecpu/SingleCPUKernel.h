//
// Created by clonker on 07.03.16.
//

#ifndef READDY2_MAIN_SINGLECPUKERNEL_H
#define READDY2_MAIN_SINGLECPUKERNEL_H

#include <Kernel.h>

namespace readdy {
    namespace kernel {
        class SingleCPUKernel : public readdy::plugin::Kernel {
        public:
            SingleCPUKernel() : Kernel("SingleCPU") { }

        };
    }
}

#endif //READDY2_MAIN_SINGLECPUKERNEL_H
