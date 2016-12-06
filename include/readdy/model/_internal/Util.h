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
 * Utility methods to deal with particle type ids: Transform a vector of types into a vector or set of type ids.
 *
 * @file Util.h
 * @brief Some utility methods for the model module.
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
