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
 * << detailed description >>
 *
 * @file ParticleTypeRegistry.h
 * @brief << brief description >>
 * @author clonker
 * @date 29.03.17
 * @copyright GNU Lesser General Public License v3.0
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
    particle_type_type typeId;

    ParticleTypeInfo(const std::string &name, scalar diffusionConstant,
                     particle_flavor flavor, Particle::type_type typeId);
};

class ParticleTypeRegistry {
public:

    using type_map = std::unordered_map<std::string, particle_type_type>;

    ParticleTypeRegistry() = default;

    ParticleTypeRegistry(const ParticleTypeRegistry &) = default;

    ParticleTypeRegistry &operator=(const ParticleTypeRegistry &) = default;

    ParticleTypeRegistry(ParticleTypeRegistry &&) = default;

    ParticleTypeRegistry &operator=(ParticleTypeRegistry &&) = default;

    ~ParticleTypeRegistry() = default;

    particle_type_type id_of(const std::string &name) const;

    particle_type_type operator()(const std::string &name) const;

    void add(const std::string &name, scalar diffusionConst, particle_flavor flavor = particleflavor::NORMAL);

    const ParticleTypeInfo &info_of(const std::string &name) const;

    const ParticleTypeInfo &info_of(Particle::type_type type) const;

    scalar diffusion_constant_of(const std::string &particleType) const;

    scalar diffusion_constant_of(particle_type_type particleType) const;

    const std::size_t &n_types() const;

    std::vector<particle_type_type> types_flat() const;

    std::string name_of(particle_type_type id) const;

    const type_map &type_mapping() const;

    void debug_output() const;

    void configure();

private:

    particle_type_type _id_of(const std::string& name) const;

    std::size_t n_types_ = 0;
    particle_type_type type_counter_ = 0;
    type_map type_mapping_ {};
    std::unordered_map<particle_type_type, ParticleTypeInfo> particle_info_ {};

};

NAMESPACE_END(model)
NAMESPACE_END(readdy)
