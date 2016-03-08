//
// Created by clonker on 07.03.16.
//

#include <iostream>
#include <Kernel.h>
#include <boost/filesystem.hpp>

namespace fs = boost::filesystem;

readdy::plugin::KernelProvider &readdy::plugin::KernelProvider::getInstance() {
    static readdy::plugin::KernelProvider instance;
    // TODO initialize kernels (load by directory) -- use boost dll
    return instance;
}

readdy::plugin::KernelProvider::KernelProvider() {
    fs::path path = fs::current_path();
    std::cout << "current path is " << path << std::endl;
}

const std::string readdy::plugin::KernelProvider::getDefaultKernelDirectory() {
    // TODO
    return "";
}

void readdy::plugin::KernelProvider::loadKernelsFromDirectory(std::string directory) {
    // TODO
}

const std::string readdy::plugin::Kernel::getName() {
    return this->name;
}

readdy::plugin::Kernel::Kernel(std::string name) {
    this->name = name;
    BOOST_LOG_TRIVIAL(trace) << "creating kernel " << name;
}
