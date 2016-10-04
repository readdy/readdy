/**
 * << detailed description >>
 *
 * @file Util.h.h
 * @brief << brief description >>
 * @author clonker
 * @date 09.08.16
 */

#ifndef READDY_MAIN_UTIL_H_H
#define READDY_MAIN_UTIL_H_H

#include <readdy/model/KernelContext.h>

namespace readdy {
namespace model {
namespace _internal {
namespace util {

inline std::set<unsigned int> transformTypes(std::vector<std::string> types, const readdy::model::KernelContext &ctx) {
    std::set<unsigned int> result;
    for (auto &&t : types) {
        result.insert(ctx.getParticleTypeID(t));
    }
    return result;
}

inline std::vector<unsigned int>
transformTypes2(std::vector<std::string> types, const readdy::model::KernelContext &ctx) {
    std::vector<unsigned int> result;
    result.reserve(types.size());
    for (auto &&t : types) {
        result.push_back(ctx.getParticleTypeID(t));
    }
    return result;
}

}
}
}
}

#endif //READDY_MAIN_UTIL_H_H
