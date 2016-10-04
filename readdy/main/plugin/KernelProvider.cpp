//
// Created by clonker on 07.03.16.
//

#include <iostream>
#include <boost/dll.hpp>
#include <boost/range/iterator_range_core.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/exception/diagnostic_information.hpp>
#include <readdy/common/Utils.h>
#include <readdy/common/make_unique.h>
#include <readdy/model/Kernel.h>
#include <readdy/plugin/KernelProvider.h>
#include <readdy/plugin/_internal/KernelPluginDecorator.h>
#include <readdy/kernel/singlecpu/SingleCPUKernel.h>

namespace fs = boost::filesystem;
namespace utl = readdy::util;
namespace readdy {
namespace plugin {

KernelProvider &KernelProvider::getInstance() {
    static KernelProvider instance;
    return instance;
}

KernelProvider::KernelProvider() {
    fs::path path = fs::current_path();
    BOOST_LOG_TRIVIAL(debug) << "current path is " << path;
    add(readdy::kernel::singlecpu::SingleCPUKernel::name, [] {
        return new readdy::kernel::singlecpu::SingleCPUKernel();
    });
}

const std::string KernelProvider::getDefaultKernelDirectory() {
    auto &&dir = getenv("READDY_PLUGIN_DIR");
    static std::string defaultDir;
    if (dir == NULL) {
        if (utl::isWindows()) {
            dir = getenv("PROGRAMFILES");
            if (dir == NULL) {
                defaultDir = "C:\\\\Program Files\\ReaDDy\\readdy\\readdy_plugins";
            } else {
                defaultDir = std::string(dir).append("\\ReaDDy\\readdy\\readdy_plugins");
            }
        } else {
            defaultDir = "/usr/local/readdy/readdy/readdy_plugins";
        }
    } else {
        defaultDir = std::string(dir);
    }
    return defaultDir;
}

void KernelProvider::loadKernelsFromDirectory(const std::string &directory) {
    BOOST_LOG_TRIVIAL(debug) << "loading kernels from directory: " << directory;
    const fs::path p(directory);
    if (fs::exists(p) && fs::is_directory(p)) {
        for (auto &&dirEntry : boost::make_iterator_range(fs::directory_iterator(p), {})) {
            try {
                if (isSharedLibrary(dirEntry.path())) {
                    add(dirEntry.path());
                }
            } catch (const std::exception &e) {
                BOOST_LOG_TRIVIAL(warning) << "Could not load " << dirEntry << " due to: " << e.what();
            } catch (const boost::exception &e) {
                BOOST_LOG_TRIVIAL(warning) << "Could not load " << dirEntry << " due to: "
                                           << boost::diagnostic_information(e);
            }
        }
    } else {
        throw std::runtime_error("file [" + p.string() + "] did not exist or was a file.");
    }
}


bool KernelProvider::isSharedLibrary(const boost::filesystem::path &path) const {
    const std::string s = path.string();
    return (s.find(".dll") != std::string::npos || s.find(".so") != std::string::npos ||
            s.find(".dylib") != std::string::npos)
           && s.find(".lib") == std::string::npos
           && s.find(".exp") == std::string::npos
           && s.find(".pdb") == std::string::npos
           && s.find(".manifest") == std::string::npos
           && s.find(".rsp") == std::string::npos
           && s.find(".obj") == std::string::npos
           && s.find(".a") == std::string::npos;
}

void KernelProvider::add(const boost::filesystem::path &sharedLib) {
    using namespace _internal;
    const auto name = loadKernelName(sharedLib);
    factory.emplace(std::make_pair(name, [sharedLib] { return new KernelPluginDecorator(sharedLib); }));
}

void KernelProvider::add(const std::string name, const std::function<readdy::model::Kernel *()> creator) {
    PluginProvider::add(name, creator);
}
}
}




