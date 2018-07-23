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
