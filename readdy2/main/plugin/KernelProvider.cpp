//
// Created by clonker on 07.03.16.
//

#include <iostream>
#include <Kernel.h>
#include <boost/range/iterator_range_core.hpp>
#include <boost/algorithm/string.hpp>
#include <Utils.h>

namespace fs = boost::filesystem;
namespace utl = readdy::utils;

readdy::plugin::KernelProvider &readdy::plugin::KernelProvider::getInstance() {
    static readdy::plugin::KernelProvider instance;
    // TODO initialize kernels (load by directory) -- use boost dll (?)
    return instance;
}

readdy::plugin::KernelProvider::KernelProvider() {
    fs::path path = fs::current_path();
    std::cout << "current path is " << path << std::endl;
}

const std::string readdy::plugin::KernelProvider::getDefaultKernelDirectory() {
    auto dir = getenv("READDY_PLUGIN_DIR");
    static std::string defaultDir;
    if(dir == NULL) {
        if(utl::isWindows()) {
            dir = getenv("PROGRAMFILES");
            if(dir == NULL) {
                defaultDir = "C:\\\\Program Files\\ReaDDy2\\lib\\plugins";
            } else {
                defaultDir = std::string(dir).append("\\ReaDDy2\\lib\\plugins");
            }
        } else {
            defaultDir = "/usr/local/readdy2/lib/plugins";
        }
    } else {
        defaultDir = std::string(dir);
    }
    return defaultDir;
}

void readdy::plugin::KernelProvider::loadKernelsFromDirectory(std::string directory) {
    const fs::path p(directory);
    if (fs::exists(p) && fs::is_directory(p)) {
        BOOST_LOG_TRIVIAL(debug) << "attempting to load plugins from directory " << p.string();
        for(auto &dirEntry : boost::make_iterator_range(fs::directory_iterator(p), {})) {
            BOOST_LOG_TRIVIAL(debug) << "... loading " << dirEntry.path().string();
        }
    } else {
        // TODO raise
        BOOST_LOG_TRIVIAL(error) << "file [" << p.string() << "] did not exist or was a file.";
    }
}

const std::string readdy::plugin::Kernel::getName() {
    return this->name;
}

readdy::plugin::Kernel::Kernel(std::string name) {
    this->name = name;
    BOOST_LOG_TRIVIAL(trace) << "creating kernel " << name;
}
