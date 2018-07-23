/********************************************************************
 * Copyright © 2018 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * Redistribution and use in source and binary forms, with or       *
 * without modification, are permitted provided that the            *
 * following conditions are met:                                    *
 *  1. Redistributions of source code must retain the above         *
 *     copyright notice, this list of conditions and the            *
 *     following disclaimer.                                        *
 *  2. Redistributions in binary form must reproduce the above      *
 *     copyright notice, this list of conditions and the following  *
 *     disclaimer in the documentation and/or other materials       *
 *     provided with the distribution.                              *
 *  3. Neither the name of the copyright holder nor the names of    *
 *     its contributors may be used to endorse or promote products  *
 *     derived from this software without specific                  *
 *     prior written permission.                                    *
 *                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND           *
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,      *
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF         *
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE         *
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR            *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,     *
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,         *
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER *
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,      *
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)    *
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF      *
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
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




