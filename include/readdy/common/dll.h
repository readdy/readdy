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


/**
 * This header file is mainly to replace boost.dll. It contains a shared_library struct that loads and unloads shared
 * libraries and is also capable of calling exported functions thereof.
 *
 * @file dll.h
 * @brief Definition of shared_library.
 * @author clonker
 * @date 14.10.16
 */

#pragma once

#include <string>
#if READDY_OSX || READDY_LINUX
#include <dlfcn.h>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(util)
NAMESPACE_BEGIN(dll)

struct shared_library {

    void* handle = nullptr;

    shared_library(const std::string& path, int mode) {
        handle = dlopen(path.c_str(), mode);
    }

    shared_library(const shared_library&) = delete;
    shared_library& operator=(const shared_library&) = delete;

    bool has_symbol(const std::string& symbol) {
        if(handle) {
            return dlsym(handle, symbol.c_str()) != nullptr;
        }
        return false;
    }

    virtual ~shared_library() {
        if(handle) {
            dlclose(handle);
            handle = 0;
        }
    }

    template<typename Signature>
    std::function<Signature> load(const std::string &symbol) {
        dlerror();
        const auto result = dlsym(handle, symbol.c_str());
        if(!result) {
            const auto error = dlerror();
            if(error) {
                throw std::logic_error("couldn't find symbol \""+symbol+"\": " + error);
            }
        }
        return reinterpret_cast<Signature*>(result);
    }

};

NAMESPACE_END(dll)
NAMESPACE_END(util)
NAMESPACE_END(readdy)

#endif
