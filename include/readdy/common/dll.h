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

    template<typename Ret, typename... Args>
    std::function<Ret(Args...)> call(const std::string& symbol, Args... args) {
        using fptr = Ret (*)(Args...);
        if(so_ptr) {
            return (fptr) (dlsym(so_ptr, symbol.c_str()));
        }
        /**
         * TODO!
         */
        /*

        if(so_ptr) {
            fptr fptr = ;
            return fptr(args...);
        }*/
        throw std::runtime_error("no library loaded");
    };

};
}
}
}
#endif //READDY_MAIN_DLL_H
