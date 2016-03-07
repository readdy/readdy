//
// Created by clonker on 07.03.16.
//

#include <Kernel.h>

std::unique_ptr<readdy::plugin::Kernel> readdy::plugin::KernelFactory::create(std::string name) {
    // TODO implement this
    return nullptr;
}

readdy::plugin::KernelFactory &readdy::plugin::KernelFactory::getInstance() {
    static readdy::plugin::KernelFactory instance;
    // TODO initialize kernels (load by directory)
    return instance;
}

readdy::plugin::KernelFactory::KernelFactory() {
    // no op
}
