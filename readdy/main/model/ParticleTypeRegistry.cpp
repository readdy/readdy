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
 * @file ParticleTypeRegistry.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 29.03.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <readdy/model/ParticleTypeRegistry.h>

namespace readdy {
namespace model {


ParticleTypeInfo::ParticleTypeInfo(const std::string &name, const double diffusionConstant, const double radius,
                                   const Particle::flavor_t flavor, const Particle::type_type typeId)
        : name(name), diffusionConstant(diffusionConstant), radius(radius), flavor(flavor), typeId(typeId) {}


const ParticleTypeInfo &ParticleTypeRegistry::getParticleTypeInfo(const std::string &name) const {
    return getParticleTypeInfo(typeMapping.at(name));
}

const ParticleTypeInfo &ParticleTypeRegistry::getParticleTypeInfo(const Particle::type_type type) const {
    return particleInfo.at(type);
}

const ParticleTypeRegistry::rdy_reverse_type_mapping ParticleTypeRegistry::generateReverseTypeMapping() const {
    const auto &typeMapping = getTypeMapping();
    rdy_reverse_type_mapping reverseTypeMapping;
    auto it = typeMapping.cbegin();
    while (it != typeMapping.cend()) {
        auto key = (*it).first;
        auto value = (*it).second;
        reverseTypeMapping.emplace(std::make_pair(value, key));
        ++it;
    }
    return reverseTypeMapping;
}

const ParticleTypeRegistry::rdy_type_mapping &ParticleTypeRegistry::getTypeMapping() const {
    return typeMapping;
}

std::string ParticleTypeRegistry::getParticleName(particle_type_type id) const {
    for (const auto &e : typeMapping) {
        if (e.second == id) return e.first;
    }
    return "";
}

std::vector<particle_type_type> ParticleTypeRegistry::getAllRegisteredParticleTypes() const {
    std::vector<particle_type_type> v;
    for (auto &&entry : typeMapping) {
        v.push_back(entry.second);
    }
    return v;
}

double ParticleTypeRegistry::getParticleRadius(const particle_type_type type) const {
    return particleInfo.at(type).radius;
}

double ParticleTypeRegistry::getParticleRadius(const std::string &type) const {
    return getParticleRadius(getParticleTypeID(type));
}

double ParticleTypeRegistry::getDiffusionConstant(const particle_type_type particleType) const {
    return particleInfo.at(particleType).diffusionConstant;
}

double ParticleTypeRegistry::getDiffusionConstant(const std::string &particleType) const {
    return getDiffusionConstant(getParticleTypeID(particleType));
}

void
ParticleTypeRegistry::registerParticleType(const std::string &name, const double diffusionConst, const double radius,
                                           const Particle::flavor_t flavor) {
    particle_type_type t_id = typeCounter++;
    typeMapping.emplace(name, t_id);
    particleInfo.emplace(std::make_pair(t_id, ParticleTypeInfo{name, diffusionConst, radius, flavor, t_id}));
}

particle_type_type ParticleTypeRegistry::getParticleTypeID(const std::string &name) const {
    return typeMapping.at(name);
}
}
}