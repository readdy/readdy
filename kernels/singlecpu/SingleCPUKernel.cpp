//
// Created by clonker on 07.03.16.
//

#include "SingleCPUKernel.h"

namespace kern = readdy::kernel::singlecpu;
kern::SingleCPUKernel::SingleCPUKernel() : readdy::plugin::Kernel("SingleCPU"){
    BOOST_LOG_TRIVIAL(debug) << "Single CPU Kernel instantiated!";
}

kern::SingleCPUKernel* kern::SingleCPUKernel::create() {
    //auto kernel = kern::SingleCPUKernel();
    return new kern::SingleCPUKernel();
}
