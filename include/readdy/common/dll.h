/**
 * This header file is mainly to replace boost.dll. It contains a shared_library struct that loads and unloads shared
 * libraries and is also capable of calling exported functions thereof.
 *
 * @file dll.h
 * @brief Definition of shared_library.
 * @author clonker
 * @date 14.10.16
 */

#ifndef READDY_MAIN_DLL_H
#define READDY_MAIN_DLL_H

#include <string>
#include <dlfcn.h>

namespace readdy {
namespace util {
namespace dll {
struct shared_library {

    void* handle = nullptr;

    shared_library() {
        // no op
    }

    shared_library(const std::string& path, int mode) {
        handle = dlopen(path.c_str(), mode);
    }

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
    };

};
}
}
}

#if READDY_LINUX || READDY_OSX
#define READDY_EXPORT_FUNCTION_AS(FunctionName, ExportName) extern "C" const void* ExportName();
#define READDY_EXPORT_FUNCTION_AS_IMPL(FunctionName, ExportName) \
     const void* ExportName() { return reinterpret_cast<const void*>(reinterpret_cast<intptr_t>(&FunctionName)); }
#else
// todo msvc
#endif // if readdy_linux || readdy_osx

#endif //READDY_MAIN_DLL_H
