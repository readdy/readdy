//
// Created by clonker on 07.03.16.
//

#include <iostream>
#include <Kernel.h>
#include <boost/filesystem.hpp>

namespace fs = boost::filesystem;

std::unique_ptr<readdy::plugin::Kernel> readdy::plugin::KernelProvider::get(std::string name) {
    // TODO implement this
    return nullptr;
}

readdy::plugin::KernelProvider &readdy::plugin::KernelProvider::getInstance() {
    static readdy::plugin::KernelProvider instance;
    // TODO initialize kernels (load by directory)
    return instance;
}

readdy::plugin::KernelProvider::KernelProvider() {
    fs::path path = fs::current_path();
    std::cout << "current path is " << path << std::endl;
}
