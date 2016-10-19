//
// Created by clonker on 07.03.16.
//

#include <iostream>
#include <readdy/common/Utils.h>
#include <readdy/common/make_unique.h>
#include <readdy/common/filesystem.h>
#include <readdy/model/Kernel.h>
#include <readdy/plugin/KernelProvider.h>
#include <readdy/plugin/_internal/KernelPluginDecorator.h>
#include <readdy/kernel/singlecpu/SingleCPUKernel.h>

namespace utl = readdy::util;
namespace fs = utl::fs;
namespace readdy {
namespace plugin {

KernelProvider &KernelProvider::getInstance() {
    static KernelProvider instance;
    return instance;
}

KernelProvider::KernelProvider() {
    const auto path = fs::current_path();
    if(!log::console()) {
        spdlog::set_sync_mode();
        auto console = spdlog::stdout_color_mt("console");
        console->set_level(spdlog::level::debug);
        console->set_pattern("[          ] [%Y-%m-%d %H:%M:%S] [%t] [%l] %v");
        log::console()->warn("initialized default console logger because there was none");
    }
    log::console()->debug("current path is {}", path);
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
    log::console()->debug("loading kernels from directory: {}", directory);
    if (fs::exists(directory) && fs::is_directory(directory)) {
        fs::dir_iterator dir(directory);
        while(dir.has_next()) {
            const auto file = dir.next();
            if(fs::is_file(file)) {
                try {
                    if (isSharedLibrary(file)) {
                        add(file);
                    }
                } catch (const std::exception &e) {
                    log::console()->warn("could not load {} due to {}", file, std::string(e.what()));
                }
            } else {
                log::console()->debug("omitted {} because it was no file.", file);
            }
        }
    } else {
        throw std::runtime_error("file [" + directory + "] did not exist or was a file.");
    }
}


bool KernelProvider::isSharedLibrary(const std::string &path) const {
    return (path.find(".dll") != std::string::npos || path.find(".so") != std::string::npos ||
            path.find(".dylib") != std::string::npos)
           && path.find(".lib") == std::string::npos
           && path.find(".exp") == std::string::npos
           && path.find(".pdb") == std::string::npos
           && path.find(".manifest") == std::string::npos
           && path.find(".rsp") == std::string::npos
           && path.find(".obj") == std::string::npos
           && path.find(".a") == std::string::npos;
}

void KernelProvider::add(const std::string &sharedLib) {
    using namespace _internal;
    const auto name = loadKernelName(sharedLib);
    log::console()->debug("Trying to load kernel with name {}", name);
    factory.emplace(std::make_pair(name, [sharedLib] { return new KernelPluginDecorator(sharedLib); }));
}

void KernelProvider::add(const std::string name, const std::function<readdy::model::Kernel *()> creator) {
    PluginProvider::add(name, creator);
}
}
}




