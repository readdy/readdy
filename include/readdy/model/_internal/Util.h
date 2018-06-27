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
 * Utility methods to deal with particle type ids: Transform a container of type strings to a container of type ids.
 *
 * @file Util.h
 * @brief Some utility methods for the model module.
 * @author clonker
 * @author chrisfroe
 * @date 09.08.16
 */

#pragma once

#include <set>
#include <sstream>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(_internal)
NAMESPACE_BEGIN(util)

template<typename Context>
inline std::set<ParticleTypeId> transformTypes(const std::vector<std::string> &types, const Context &ctx) {
    std::set<ParticleTypeId> result;
    for (const auto &t : types) {
        result.insert(ctx.particleTypes().idOf(t));
    }
    return result;
}

template<typename Context>
inline std::vector<ParticleTypeId>
transformTypes2(const std::vector<std::string> &types, const Context &ctx) {
    std::vector<ParticleTypeId> result;
    result.reserve(types.size());
    for (auto &t : types) {
        result.push_back(ctx.particleTypes().idOf(t));
    }
    return result;
}

inline std::unordered_map<Particle::type_type, Particle::type_type>
transformTypesMap(const std::unordered_map<std::string, std::string> &stringMap, const ParticleTypeRegistry &types) {
    std::unordered_map<Particle::type_type, Particle::type_type> result;
    for (const auto &pair : stringMap) {
        result.emplace(types(pair.first), types(pair.second));
    }
    return result;
}

template<typename T>
std::string to_string(const T& ref) {
    std::stringstream ss;
    ss << ref;
    return ss.str();
}

template<typename T>
std::string to_string(const T* const ptr) {
    std::stringstream ss;
    ss << *ptr;
    return ss.str();
}

NAMESPACE_END(util)
NAMESPACE_END(_internal)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
