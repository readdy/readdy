//
// Created by clonker on 07.03.16.
//

#include <iostream>
#include <Kernel.h>
#include <boost/range/iterator_range_core.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/dll.hpp>
#include <Utils.h>

namespace fs = boost::filesystem;
namespace utl = readdy::utils;
namespace plug = readdy::plugin;

plug::KernelProvider &plug::KernelProvider::getInstance() {
    static plug::KernelProvider instance;
    // TODO initialize kernels (load by directory) -- use boost dll (?)
    return instance;
}

plug::KernelProvider::KernelProvider() {
    fs::path path = fs::current_path();
    std::cout << "current path is " << path << std::endl;
}

const std::string plug::KernelProvider::getDefaultKernelDirectory() {
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

void plug::KernelProvider::loadKernelsFromDirectory(std::string directory) {
    const fs::path p(directory);
    typedef std::shared_ptr<plug::Kernel> (kernel_t)();
    if (fs::exists(p) && fs::is_directory(p)) {
        BOOST_LOG_TRIVIAL(debug) << "attempting to load plugins from directory " << p.string();
        BOOST_LOG_TRIVIAL(debug) << "current path: " << fs::current_path().string();
        for(auto dirEntry : boost::make_iterator_range(fs::directory_iterator(p), {})) {
            BOOST_LOG_TRIVIAL(debug) << "... loading " << dirEntry.path().string();
            boost::function<kernel_t> factory;
            factory = boost::dll::import_alias<kernel_t>(
                    dirEntry.path(),
                    "create_kernel",
                    boost::dll::load_mode::rtld_lazy | boost::dll::load_mode::rtld_global
            );
            auto kernel = factory();
            plug::KernelProvider::getInstance().add(kernel);
        }
    } else {
        // TODO raise
        BOOST_LOG_TRIVIAL(error) << "file [" << p.string() << "] did not exist or was a file.";
    }
}

const std::string plug::Kernel::getName() {
    return this->name;
}

plug::Kernel::Kernel(std::string name) {
    this->name = name;
    BOOST_LOG_TRIVIAL(trace) << "creating kernel " << name;
}
