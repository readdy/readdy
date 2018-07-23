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
 * This header file contains:
 *     * the `readdy::model::particle_flavor` type which takes values as defined in the `readdy::model::particleflavor`
 *       namespace
 *     * the definition of ParticleTypeInfo, backing structure for all particle types
 *     * the definition of ParticleTypeRegistry, an object managing all particle types for a reaction diffusion system
 *
 * @file ParticleTypeRegistry.h
 * @brief Header file containing definitions for particle flavors, particle type info and the particle type registry.
 * @author clonker
 * @date 29.03.17
 * @copyright GPL-3
 */

#pragma once

#include <readdy/common/common.h>

#include "Particle.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)

using particle_flavor = std::uint8_t;
NAMESPACE_BEGIN(particleflavor)
static constexpr particle_flavor NORMAL = 0;
static constexpr particle_flavor TOPOLOGY = 1;
static constexpr particle_flavor MEMBRANE = 2;

inline static std::string particle_flavor_to_str(particle_flavor flavor) {
    switch(flavor) {
        case model::particleflavor::NORMAL: return "NORMAL";
        case model::particleflavor::TOPOLOGY: return "TOPOLOGY";
        case model::particleflavor::MEMBRANE: return "MEMBRANE";
        default: return "UNKNOWN";
    }
}
NAMESPACE_END(particleflavor)

struct ParticleTypeInfo {
    std::string name;
    scalar diffusionConstant;
    particle_flavor flavor;
    ParticleTypeId typeId;

    ParticleTypeInfo(const std::string &name, scalar diffusionConstant,
                     particle_flavor flavor, Particle::type_type typeId);
};

class ParticleTypeRegistry {
public:

    using type_map = std::unordered_map<std::string, ParticleTypeId>;

    ParticleTypeRegistry() = default;

    ParticleTypeRegistry(const ParticleTypeRegistry &) = default;

    ParticleTypeRegistry &operator=(const ParticleTypeRegistry &) = default;

    ParticleTypeRegistry(ParticleTypeRegistry &&) = default;

    ParticleTypeRegistry &operator=(ParticleTypeRegistry &&) = default;

    ~ParticleTypeRegistry() = default;

    ParticleTypeId idOf(const std::string &name) const {
        return _idOf(name);
    }

    ParticleTypeId operator()(const std::string &name) const {
        return idOf(name);
    }

    void add(const std::string &name, scalar diffusionConst, particle_flavor flavor = particleflavor::NORMAL);

    void addTopologyType(const std::string &name, scalar diffusionConst) {
        add(name, diffusionConst, particleflavor::TOPOLOGY);
    }

    const ParticleTypeInfo &infoOf(const std::string &name) const {
        return infoOf(_idOf(name));
    }

    const ParticleTypeInfo &infoOf(Particle::type_type type) const {
        return particle_info_.at(type);
    }

    scalar diffusionConstantOf(const std::string &particleType) const {
        return diffusionConstantOf(idOf(particleType));
    }

    scalar& diffusionConstantOf(const std::string &particleType) {
        return diffusionConstantOf(idOf(particleType));
    }

    scalar diffusionConstantOf(ParticleTypeId particleType) const {
        return particle_info_.at(particleType).diffusionConstant;
    }

    scalar& diffusionConstantOf(ParticleTypeId particleType) {
        return particle_info_.at(particleType).diffusionConstant;
    }

    const std::size_t &nTypes() const {
        return n_types_;
    }

    std::vector<ParticleTypeId> typesFlat() const {
        std::vector<ParticleTypeId> v;
        std::transform(std::begin(type_mapping_), std::end(type_mapping_), std::back_inserter(v), [](const auto &e) {
            return e.second;
        });
        return v;
    }

    std::string nameOf(ParticleTypeId id) const {
        auto it = std::find_if(std::begin(type_mapping_), std::end(type_mapping_), [id](const auto &e) {
            return e.second == id;
        });
        return it == std::end(type_mapping_) ? "" : it->first;
    }

    const type_map &typeMapping() const {
        return type_mapping_;
    }

    std::string describe() const;

private:

    ParticleTypeId _idOf(const std::string &name) const {
        auto it = type_mapping_.find(name);
        if (it == type_mapping_.end()) {
            throw std::invalid_argument(
                    fmt::format("Could not find type \"{}\", did you forget to register it before accessing it?", name)
            );
        }
        return it->second;
    }

    std::size_t n_types_ = 0;
    ParticleTypeId type_counter_ = 0;
    type_map type_mapping_ {};
    std::unordered_map<ParticleTypeId, ParticleTypeInfo> particle_info_ {};

};

NAMESPACE_END(model)
NAMESPACE_END(readdy)
