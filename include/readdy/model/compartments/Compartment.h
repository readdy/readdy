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
 * << detailed description >>
 *
 * @file Compartment.h
 * @brief << brief description >>
 * @author chrisfroe
 * @date 12.01.17
 * @copyright BSD-3
 */
#pragma once

#include <unordered_map>
#include <utility>

#include <readdy/model/Particle.h>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(compartments)

class Compartment {
public:
    using id_type = short;
    using label_conversion_map = std::unordered_map<std::string, std::string>;
    using conversion_map = std::unordered_map<ParticleTypeId, ParticleTypeId>;

    Compartment(conversion_map conversions, std::string typeName,
                std::string uniqueName)
            : typeName(std::move(typeName)), uniqueName(std::move(uniqueName)),
              conversions(std::move(conversions)), _id(counter++) {}

    virtual const bool isContained(const Vec3 &position) const = 0;

    const conversion_map &getConversions() const {
        return conversions;
    }

    const std::string &getTypeName() const {
        return typeName;
    }

    const std::string &getUniqueName() const {
        return uniqueName;
    }

    const id_type getId() const {
        return _id;
    }

protected:
    static id_type counter;

    std::string typeName;
    std::string uniqueName;
    id_type _id;
    conversion_map conversions;
};

NAMESPACE_END(compartments)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
