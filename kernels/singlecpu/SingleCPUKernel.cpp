//
// Created by clonker on 07.03.16.
//

#include "SingleCPUKernel.h"

namespace kern = readdy::kernel::singlecpu;
kern::SingleCPUKernel::SingleCPUKernel() : readdy::plugin::Kernel("SingleCPU"){
    BOOST_LOG_TRIVIAL(debug) << "Single CPU Kernel instantiated!";
}

std::shared_ptr<kern::SingleCPUKernel> kern::SingleCPUKernel::create() {
    return std::shared_ptr<kern::SingleCPUKernel>(new kern::SingleCPUKernel());
}
