/**
 * << detailed description >>
 *
 * @file dll.h
 * @brief << brief description >>
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

    void* so_ptr = nullptr;

    shared_library() {
        // no op
    }

    shared_library(const std::string& path, int mode) {
        so_ptr = dlopen(path.c_str(), mode);
    }

    bool has_symbol(const std::string& symbol) {
        if(so_ptr) {
            return dlsym(so_ptr, symbol.c_str()) != nullptr;
        }
        return false;
    }

    virtual ~shared_library() {
        if(so_ptr) {
            dlclose(so_ptr);
        }
    }

    template<typename Signature>
    std::function<Signature> load(const std::string &symbol) {
        dlerror();
        const auto result = dlsym(so_ptr, symbol.c_str());
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
#endif //READDY_MAIN_DLL_H
