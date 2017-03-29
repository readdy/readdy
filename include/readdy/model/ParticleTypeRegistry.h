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

struct ParticleTypeInfo {
    std::string name;
    double diffusionConstant;
    double radius;
    readdy::model::Particle::flavor_t flavor;
    readdy::model::Particle::type_type typeId;

    ParticleTypeInfo(const std::string &name, const double diffusionConstant, const double radius,
                     const Particle::flavor_t flavor, const Particle::type_type typeId);
};

class ParticleTypeRegistry {
public:

    using rdy_type_mapping = std::unordered_map<std::string, particle_type_type>;
    using rdy_reverse_type_mapping = std::unordered_map<particle_type_type, std::string>;

    ParticleTypeRegistry() = default;
    ParticleTypeRegistry(const ParticleTypeRegistry&) = delete;
    ParticleTypeRegistry& operator=(const ParticleTypeRegistry&) = delete;
    ParticleTypeRegistry(ParticleTypeRegistry&&) = delete;
    ParticleTypeRegistry& operator=(ParticleTypeRegistry&&) = delete;

    particle_type_type getParticleTypeID(const std::string &name) const;

    void registerParticleType(const std::string &name, const double diffusionConst, const double radius,
                              const readdy::model::Particle::flavor_t flavor = readdy::model::Particle::FLAVOR_NORMAL);

    const ParticleTypeInfo &getParticleTypeInfo(const std::string &name) const;

    const ParticleTypeInfo &getParticleTypeInfo(const Particle::type_type type) const;

    double getDiffusionConstant(const std::string &particleType) const;

    double getDiffusionConstant(const particle_type_type particleType) const;

    double getParticleRadius(const std::string &type) const;

    double getParticleRadius(const particle_type_type type) const;

    std::vector<particle_type_type> getAllRegisteredParticleTypes() const;

    std::string getParticleName(particle_type_type id) const;

    const rdy_type_mapping &getTypeMapping() const;

    /**
     * Generate a map from particle_type_t to string. As there is no book-keeping of this reversed
     * structure, it is generated in place and then returned.
     */
    const rdy_reverse_type_mapping generateReverseTypeMapping() const;

private:
    particle_type_type typeCounter = 0;
    rdy_type_mapping typeMapping;
    std::unordered_map<particle_type_type, ParticleTypeInfo> particleInfo;

};

NAMESPACE_END(model)
NAMESPACE_END(readdy)
