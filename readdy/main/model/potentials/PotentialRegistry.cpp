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


#include <readdy/model/potentials/PotentialRegistry.h>
#include <readdy/common/Utils.h>
#include <readdy/model/potentials/PotentialsOrder1.h>
#include <readdy/common/string.h>

/**
 * << detailed description >>
 *
 * @file PotentialRegistry.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 29.03.17
 * @copyright GNU Lesser General Public License v3.0
 */

namespace readdy {
namespace model {
namespace potentials {
PotentialRegistry::PotentialRegistry(particle_type_registry_ref typeRegistry) : _types(typeRegistry) {}

PotentialRegistry::id_type PotentialRegistry::addUserDefined(potentials::PotentialOrder2 *potential) {
    auto id = potential->getId();
    auto type1Id = potential->particleType1();
    auto type2Id = potential->particleType2();
    auto pp = std::tie(type1Id, type2Id);
    if (potentialO2RegistryExternal.find(pp) == potentialO2RegistryExternal.end()) {
        potentialO2RegistryExternal.emplace(pp, pot_ptr_vec2_external());
    }
    potentialO2RegistryExternal[pp].push_back(potential);
    return id;
}

PotentialRegistry::id_type PotentialRegistry::addUserDefined(potentials::PotentialOrder1 *potential) {
    auto typeId = potential->particleType();
    if (potentialO1RegistryExternal.find(typeId) == potentialO1RegistryExternal.end()) {
        potentialO1RegistryExternal.emplace(std::make_pair(typeId, pot_ptr_vec1_external()));
    }
    potentialO1RegistryExternal[typeId].push_back(potential);
    return potential->getId();
}

void PotentialRegistry::remove(const Potential::id_type handle) {
    for (auto &entry : potentialO1RegistryInternal) {
        entry.second.erase(std::remove_if(entry.second.begin(), entry.second.end(),
                                          [=](const std::shared_ptr<potentials::PotentialOrder1> &p) -> bool {
                                              return handle == p->getId();
                                          }
        ), entry.second.end());
    }
    for (auto &entry : potentialO2RegistryInternal) {
        entry.second.erase(std::remove_if(entry.second.begin(), entry.second.end(),
                                          [=](const std::shared_ptr<potentials::PotentialOrder2> &p) -> bool {
                                              return handle == p->getId();
                                          }
        ), entry.second.end());
    }
    for (auto &entry : potentialO1RegistryExternal) {
        entry.second.erase(std::remove_if(entry.second.begin(), entry.second.end(),
                                          [=](potentials::PotentialOrder1 *p) -> bool {
                                              return handle == p->getId();
                                          }
        ), entry.second.end());
    }
    for (auto &entry : potentialO2RegistryExternal) {
        entry.second.erase(std::remove_if(entry.second.begin(), entry.second.end(),
                                          [=](potentials::PotentialOrder2 *p) -> bool {
                                              return handle == p->getId();
                                          }
        ), entry.second.end());
    }
}

const PotentialRegistry::potentials_o1 &PotentialRegistry::potentialsOf(const particle_type_type type) const {
    auto it = potentialO1Registry.find(type);
    if (it != potentialO1Registry.end()) {
        return it->second;
    }
    return defaultPotentialsO1;
}

const PotentialRegistry::potential_o1_registry &PotentialRegistry::potentialsOrder1() const {
    return potentialO1Registry;
}

const std::vector<PotentialOrder2 *> &
PotentialRegistry::potentialsOf(const particle_type_type t1, const particle_type_type t2) const {
    auto it = potentialO2Registry.find(std::tie(t1, t2));
    if (it != potentialO2Registry.end()) {
        return it->second;
    }
    return defaultPotentialsO2;
}

const PotentialRegistry::potential_o2_registry &PotentialRegistry::potentialsOrder2() const {
    return potentialO2Registry;
}

const PotentialRegistry::potentials_o1 &PotentialRegistry::potentialsOf(const std::string &type) const {
    return potentialsOf(_types.get().idOf(type));
}

const PotentialRegistry::potentials_o2 &
PotentialRegistry::potentialsOf(const std::string &t1, const std::string &t2) const {
    return potentialsOf(_types.get().idOf(t1), _types.get().idOf(t2));
}

void PotentialRegistry::configure() {
    namespace coll = readdy::util::collections;
    using pair = util::particle_type_pair;
    using pot1 = potentials::PotentialOrder1;
    using pot1_ptr = std::shared_ptr<potentials::PotentialOrder1>;
    using pot2_ptr = std::shared_ptr<potentials::PotentialOrder2>;
    using pot2 = potentials::PotentialOrder2;
    potentialO1Registry.clear();
    potentialO2Registry.clear();

    coll::for_each_value(potentialO1RegistryInternal, [&](const particle_type_type type, const pot1_ptr &ptr) {
        (potentialO1Registry)[type].push_back(ptr.get());
    });
    coll::for_each_value(potentialO2RegistryInternal, [&](const pair &type, const pot2_ptr &ptr) {
        (potentialO2Registry)[type].push_back(ptr.get());
    });
    coll::for_each_value(potentialO1RegistryExternal, [&](const particle_type_type type, pot1 *ptr) {
        (potentialO1Registry)[type].push_back(ptr);
    });
    coll::for_each_value(potentialO2RegistryExternal, [&](const pair &type, pot2 *ptr) {
        (potentialO2Registry)[type].push_back(ptr);
    });
}

std::string PotentialRegistry::describe() const {
    namespace rus = readdy::util::str;
    auto find_pot_name = [this](particle_type_type type) -> const std::string {
        for (auto &&t : _types.get().typeMapping()) {
            if (t.second == type) return t.first;
        }
        return "";
    };
    std::string description;
    if (!potentialsOrder1().empty()) {
        description += fmt::format(" - potentials of order 1:{}", rus::newline);
        for (const auto &types : potentialsOrder1()) {
            description += fmt::format("     * for type {}{}", find_pot_name(types.first), rus::newline);
            for (auto pot : types.second) {
                description += fmt::format("         * {}{}", pot->describe(), rus::newline);
            }
        }
    }
    if (!potentialsOrder2().empty()) {
        description += fmt::format(" - potentials of order 2:{}", rus::newline);
        for (const auto &types : potentialsOrder2()) {
            description += fmt::format("     * for types {} and {}{}", find_pot_name(std::get<0>(types.first)),
                                       find_pot_name(std::get<1>(types.first)), rus::newline);
            for (auto pot : types.second) {
                description += fmt::format("         * {}{}", pot->describe(), rus::newline);
            }
        }
    }
    return description;
}

Potential::id_type PotentialRegistry::addBox(const std::string &particleType, scalar forceConstant, const Vec3 &origin,
                                             const Vec3 &extent) {
    return addBox(_types.get()(particleType), forceConstant, origin, extent);
}

Potential::id_type PotentialRegistry::addBox(particle_type_type particleType, scalar forceConstant, const Vec3 &origin,
                                             const Vec3 &extent) {
    auto &pots = potentialO1RegistryInternal[particleType];
    pots.emplace_back(std::make_shared<Box>(particleType, forceConstant, origin, extent));
    return pots.back()->getId();
}

PotentialRegistry::id_type
PotentialRegistry::addHarmonicRepulsion(const std::string &type1, const std::string &type2, scalar forceConstant, scalar interactionDistance) {
    return addHarmonicRepulsion(_types.get()(type1), _types.get()(type2), forceConstant, interactionDistance);
}

PotentialRegistry::id_type
PotentialRegistry::addHarmonicRepulsion(particle_type_type type1, particle_type_type type2, scalar forceConstant, scalar interactionDistance) {
    auto &pots = potentialO2RegistryInternal[std::tie(type1, type2)];
    pots.emplace_back(std::make_shared<HarmonicRepulsion>(type1, type2, forceConstant, interactionDistance));
    return pots.back()->getId();
}

PotentialRegistry::id_type
PotentialRegistry::addWeakInteractionPiecewiseHarmonic(const std::string &type1, const std::string &type2,
                                                       scalar forceConstant,
                                                       const WeakInteractionPiecewiseHarmonic::Configuration &config) {
    return addWeakInteractionPiecewiseHarmonic(_types.get()(type1), _types.get()(type2), forceConstant, config);
}

PotentialRegistry::id_type
PotentialRegistry::addWeakInteractionPiecewiseHarmonic(particle_type_type type1, particle_type_type type2,
                                                       scalar forceConstant,
                                                       const WeakInteractionPiecewiseHarmonic::Configuration &config) {
    auto &pots = potentialO2RegistryInternal[std::tie(type1, type2)];
    pots.emplace_back(std::make_shared<WeakInteractionPiecewiseHarmonic>(type1, type2, forceConstant, config));
    return pots.back()->getId();
}

PotentialRegistry::id_type
PotentialRegistry::addWeakInteractionPiecewiseHarmonic(particle_type_type type1, particle_type_type type2,
                                                       scalar forceConstant, scalar desiredDist, scalar depth,
                                                       scalar cutoff) {
    WeakInteractionPiecewiseHarmonic::Configuration conf{desiredDist, depth, cutoff};
    return addWeakInteractionPiecewiseHarmonic(type1, type2, forceConstant, conf);
}

PotentialRegistry::id_type
PotentialRegistry::addWeakInteractionPiecewiseHarmonic(const std::string &type1, const std::string &type2,
                                                       scalar forceConstant, scalar desiredDist, scalar depth,
                                                       scalar cutoff) {
    return addWeakInteractionPiecewiseHarmonic(_types.get()(type1), _types.get()(type2), forceConstant, desiredDist,
                                               depth, cutoff);
}

PotentialRegistry::id_type
PotentialRegistry::addLennardJones(const std::string &type1, const std::string &type2, unsigned int m, unsigned int n,
                                   scalar cutoff, bool shift, scalar epsilon, scalar sigma) {
    return addLennardJones(_types.get()(type1), _types.get()(type2), m, n, cutoff, shift, epsilon, sigma);
}

PotentialRegistry::id_type
PotentialRegistry::addLennardJones(particle_type_type type1, particle_type_type type2, unsigned int m, unsigned int n,
                                   scalar cutoff, bool shift, scalar epsilon, scalar sigma) {
    auto &pots = potentialO2RegistryInternal[std::tie(type1, type2)];
    pots.emplace_back(std::make_shared<LennardJones>(type1, type2, m, n, cutoff, shift, epsilon, sigma));
    return pots.back()->getId();
}

PotentialRegistry::id_type
PotentialRegistry::addScreenedElectrostatics(const std::string &particleType1, const std::string &particleType2,
                                             scalar electrostaticStrength, scalar inverseScreeningDepth,
                                             scalar repulsionStrength, scalar repulsionDistance, unsigned int exponent,
                                             scalar cutoff) {
    return addScreenedElectrostatics(_types.get()(particleType1), _types.get()(particleType2), electrostaticStrength,
                                     inverseScreeningDepth, repulsionStrength, repulsionDistance, exponent, cutoff);
}

PotentialRegistry::id_type
PotentialRegistry::addScreenedElectrostatics(particle_type_type particleType1, particle_type_type particleType2,
                                             scalar electrostaticStrength, scalar inverseScreeningDepth,
                                             scalar repulsionStrength, scalar repulsionDistance, unsigned int exponent,
                                             scalar cutoff) {
    auto &pots = potentialO2RegistryInternal[std::tie(particleType1, particleType2)];
    pots.emplace_back(std::make_shared<ScreenedElectrostatics>(particleType1, particleType2, electrostaticStrength,
                                                               inverseScreeningDepth, repulsionStrength,
                                                               repulsionDistance, exponent, cutoff));
    return pots.back()->getId();
}

PotentialRegistry::id_type
PotentialRegistry::addSphereOut(const std::string &particleType, scalar forceConstant, const Vec3 &origin,
                                scalar radius) {
    return addSphereOut(_types.get()(particleType), forceConstant, origin, radius);
}

PotentialRegistry::id_type
PotentialRegistry::addSphereOut(particle_type_type particleType, scalar forceConstant, const Vec3 &origin,
                                scalar radius) {
    auto &pots = potentialO1RegistryInternal[particleType];
    pots.emplace_back(std::make_shared<SphereOut>(particleType, forceConstant, origin, radius));
    return pots.back()->getId();
}

PotentialRegistry::id_type
PotentialRegistry::addSphereIn(const std::string &particleType, scalar forceConstant, const Vec3 &origin,
                               scalar radius) {
    return addSphereIn(_types.get()(particleType), forceConstant, origin, radius);
}

PotentialRegistry::id_type
PotentialRegistry::addSphereIn(particle_type_type particleType, scalar forceConstant, const Vec3 &origin,
                               scalar radius) {
    auto &pots = potentialO1RegistryInternal[particleType];
    pots.emplace_back(std::make_shared<SphereIn>(particleType, forceConstant, origin, radius));
    return pots.back()->getId();
}

PotentialRegistry::id_type
PotentialRegistry::addSphericalBarrier(const std::string &particleType, scalar height, scalar width, const Vec3 &origin,
                                       scalar radius) {
    return addSphericalBarrier(_types.get()(particleType), height, width, origin, radius);
}

PotentialRegistry::id_type
PotentialRegistry::addSphericalBarrier(particle_type_type particleType, scalar height, scalar width, const Vec3 &origin,
                                       scalar radius) {
    auto &pots = potentialO1RegistryInternal[particleType];
    pots.emplace_back(std::make_shared<SphericalBarrier>(particleType, height, width, origin, radius));
    return pots.back()->getId();
}


}
}
}

