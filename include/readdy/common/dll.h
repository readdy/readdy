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
