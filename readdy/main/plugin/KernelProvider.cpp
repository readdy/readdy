/********************************************************************
 * Copyright © 2016 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * This file is part of ReaDDy.                                     *
 *                                                                  *
 * ReaDDy is free software: you can redistribute it and/or modify   *
 * it under the terms of the GNU Lesser General Public License as   *
 * published by the Free Software Foundation, either version 3 of   *
 * the License, or (at your option) any later version.              *
 *                                                                  *
 * This program is distributed in the hope that it will be useful,  *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of   *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    *
 * GNU Lesser General Public License for more details.              *
 *                                                                  *
 * You should have received a copy of the GNU Lesser General        *
 * Public License along with this program. If not, see              *
 * <http://www.gnu.org/licenses/>.                                  *
 ********************************************************************/


//
// Created by clonker on 07.03.16.
//

#include <iostream>
#include <memory>

#include <readdy/common/Utils.h>
#include <readdy/common/filesystem.h>
#include <readdy/model/Kernel.h>
#include <readdy/plugin/KernelProvider.h>
#include <readdy/kernel/singlecpu/SCPUKernel.h>
#include <readdy/common/dll.h>

namespace utl = readdy::util;
namespace fs = utl::fs;
namespace dll = utl::dll;
namespace readdy {
namespace plugin {

struct KernelProvider::Impl {
    std::unordered_map<std::string, std::weak_ptr<readdy::util::dll::shared_library>> libs {};
    std::unordered_map<std::string, std::function<KernelProvider::kernel_ptr()>> factory {};
};

KernelProvider &KernelProvider::getInstance() {
    static KernelProvider instance;
    return instance;
}

KernelProvider::KernelProvider() : pimpl(std::make_unique<Impl>()){
    const auto path = fs::current_path();
    log::debug("current path is {}", path);
    add(readdy::kernel::scpu::SCPUKernel::name, [] {
        return new readdy::kernel::scpu::SCPUKernel();
    });
}

const std::string KernelProvider::getDefaultKernelDirectory() {
    auto &&dir = getenv("READDY_PLUGIN_DIR");
    static std::string defaultDir;
    if (dir == nullptr) {
        if (utl::isWindows()) {
            dir = getenv("PROGRAMFILES");
            if (dir == nullptr) {
                defaultDir = R"(C:\\Program Files\ReaDDy\readdy\readdy_plugins)";
            } else {
                defaultDir = std::string(dir).append(R"(\ReaDDy\readdy\readdy_plugins)");
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
    log::debug("loading kernels from directory: {}", directory);
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
                    log::warn("could not load {} due to {}", file, std::string(e.what()));
                }
            } else {
                log::debug("omitted {} because it was no file.", file);
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

const std::string readdy::plugin::KernelProvider::loadKernelName(const std::string &sharedLib) {
    auto it = pimpl->libs.find(sharedLib);
    if(it != pimpl->libs.end() && !it->second.expired()) {
        return it->second.lock()->load<const char *()>("name")();
    }

    dll::shared_library lib {sharedLib, RTLD_LAZY | RTLD_GLOBAL};
    if(!lib.has_symbol("name")) {
        throw std::invalid_argument("library " + sharedLib + " had no \"name\" symbol");
    } else {
        if(!lib.has_symbol("createKernel")) {
            throw std::invalid_argument("library " + sharedLib + " had no \"createKernel\" symbol");
        }
        return lib.load<const char*()>("name")();
    }
}

void KernelProvider::add(const std::string &sharedLib) {
    const auto name = loadKernelName(sharedLib);
    pimpl->factory.emplace(std::make_pair(name, [this, sharedLib, name] { // load the library
        log::debug("Trying to load kernel with name {}", name);
        auto it = pimpl->libs.find(sharedLib);
        if(it != pimpl->libs.end() && !it->second.expired()) {
            auto libPtr = it->second.lock();
            // load the kernel
            readdy::model::Kernel* kernel = libPtr->load<readdy::model::Kernel*()>("createKernel")();
            log::debug("loaded kernel with name {}", kernel->name());
            return kernel_ptr(kernel, KernelDeleter{libPtr});
        }

        auto lib = std::make_shared<readdy::util::dll::shared_library>(sharedLib, RTLD_LAZY | RTLD_GLOBAL);
        pimpl->libs[sharedLib] = lib;
        // load the kernel
        readdy::model::Kernel* kernel = lib->load<readdy::model::Kernel*()>("createKernel")();
        log::debug("loaded kernel with name {}", kernel->name());
        return kernel_ptr(kernel, KernelDeleter{lib});
    }));
}

void KernelProvider::add(const std::string &name, const std::function<readdy::model::Kernel *()> &creator) {
    pimpl->factory.emplace(std::make_pair(name, [creator]() {
        return kernel_ptr(creator(), KernelDeleter{});
    }));
}

KernelProvider::kernel_ptr KernelProvider::create(const std::string &name) const {
    auto it = pimpl->factory.find(name);
    if (it != pimpl->factory.end()) {
        return it->second();
    }
    throw std::invalid_argument("Could not load plugin with name \"" + name + "\"");
}

std::vector<std::string> KernelProvider::availableKernels() const {
    std::vector<std::string> result;
    result.reserve(pimpl->factory.size());
    for(const auto &entry : pimpl->factory) {
        result.push_back(entry.first);
    }
    return result;
}

KernelProvider::~KernelProvider() = default;

KernelDeleter::KernelDeleter()  : ptr(nullptr) {}

KernelDeleter::KernelDeleter(const std::shared_ptr<readdy::util::dll::shared_library> &libPtr) : ptr(libPtr) {}

void KernelDeleter::operator()(readdy::model::Kernel *k) {
    delete k;
    ptr.reset();
}

}
}




