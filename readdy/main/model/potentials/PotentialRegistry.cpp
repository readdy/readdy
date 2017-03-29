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
PotentialRegistry::PotentialRegistry(particle_type_registry_ref typeRegistry) : typeRegistry(typeRegistry) {}

const short PotentialRegistry::registerExternalPotential(potentials::PotentialOrder2 *potential) {
    const auto id = potential->getId();
    auto type1Id = typeRegistry.getParticleTypeID(potential->particleType1);
    auto type2Id = typeRegistry.getParticleTypeID(potential->particleType2);
    auto pp = std::tie(type1Id, type2Id);
    if (potentialO2RegistryExternal.find(pp) == potentialO2RegistryExternal.end()) {
        potentialO2RegistryExternal.emplace(pp, pot_ptr_vec2_external());
    }
    potentialO2RegistryExternal[pp].push_back(potential);
    return id;
}

const short PotentialRegistry::registerExternalPotential(potentials::PotentialOrder1 *potential) {
    auto typeId = typeRegistry.getParticleTypeID(potential->particleType);
    if (potentialO1RegistryExternal.find(typeId) == potentialO1RegistryExternal.end()) {
        potentialO1RegistryExternal.emplace(std::make_pair(typeId, pot_ptr_vec1_external()));
    }
    potentialO1RegistryExternal[typeId].push_back(potential);
    return potential->getId();
}

void PotentialRegistry::deregisterPotential(const short handle) {
    for (auto it = potentialO1RegistryInternal.begin(); it != potentialO1RegistryInternal.end(); ++it) {
        it->second.erase(std::remove_if(it->second.begin(), it->second.end(),
                                        [=](const std::unique_ptr<potentials::PotentialOrder1> &p) -> bool {
                                            return handle == p->getId();
                                        }
        ), it->second.end());
    }
    for (auto it = potentialO2RegistryInternal.begin(); it != potentialO2RegistryInternal.end(); ++it) {
        it->second.erase(std::remove_if(it->second.begin(), it->second.end(),
                                        [=](const std::unique_ptr<potentials::PotentialOrder2> &p) -> bool {
                                            return handle == p->getId();
                                        }
        ), it->second.end());
    }
    for (auto it = potentialO1RegistryExternal.begin(); it != potentialO1RegistryExternal.end(); ++it) {
        it->second.erase(std::remove_if(it->second.begin(), it->second.end(),
                                        [=](potentials::PotentialOrder1 *p) -> bool {
                                            return handle == p->getId();
                                        }
        ), it->second.end());
    }
    for (auto it = potentialO2RegistryExternal.begin(); it != potentialO2RegistryExternal.end(); ++it) {
        it->second.erase(std::remove_if(it->second.begin(), it->second.end(),
                                        [=](potentials::PotentialOrder2 *p) -> bool {
                                            return handle == p->getId();
                                        }
        ), it->second.end());
    }
}

const PotentialRegistry::potentials_o1& PotentialRegistry::getOrder1Potentials(const particle_type_type type) const {
    auto it = potentialO1Registry.find(type);
    if(it != potentialO1Registry.end()) {
        return it->second;
    }
    return defaultPotentialsO1;
}

const PotentialRegistry::rdy_pot_1_registry& PotentialRegistry::getAllOrder1Potentials() const {
    return potentialO1Registry;
}

std::vector<particle_type_type> PotentialRegistry::getAllOrder1RegisteredPotentialTypes() const {
    std::vector<particle_type_type> result{};
    result.reserve(potentialO1Registry.size());
    for (auto it = potentialO1Registry.begin(); it != potentialO1Registry.end(); ++it) {
        result.push_back(it->first);
    }
    return result;
}

std::vector<util::particle_type_pair>
PotentialRegistry::getAllOrder2RegisteredPotentialTypes() const {
    std::vector<util::particle_type_pair> result{};
    result.reserve(potentialO2Registry.size());
    for (auto it = potentialO2Registry.begin(); it != potentialO2Registry.end(); ++it) {
        result.push_back(it->first);
    }
    return result;
}

const std::vector<PotentialOrder2 *> &
PotentialRegistry::getOrder2Potentials(const particle_type_type t1, const particle_type_type t2) const {
    auto it = potentialO2Registry.find(std::tie(t1, t2));
    if(it != potentialO2Registry.end()) {
        return it->second;
    }
    return defaultPotentialsO2;
}

const PotentialRegistry::rdy_pot_2_registry &PotentialRegistry::getAllOrder2Potentials() const {
    return potentialO2Registry;
}

const PotentialRegistry::potentials_o2 PotentialRegistry::getVectorAllOrder2Potentials() const {
    std::vector<potentials::PotentialOrder2 *> result;
    for (auto &&e : getAllOrder2RegisteredPotentialTypes()) {
        for (auto &&p : getOrder2Potentials(std::get<0>(e), std::get<1>(e))) {
            result.push_back(p);
        }
    }
    return result;
}

const PotentialRegistry::potentials_o1& PotentialRegistry::getOrder1Potentials(const std::string &type) const {
    return getOrder1Potentials(typeRegistry.getParticleTypeID(type));
}

const PotentialRegistry::potentials_o2 &
PotentialRegistry::getOrder2Potentials(const std::string &t1, const std::string &t2) const {
    return getOrder2Potentials(typeRegistry.getParticleTypeID(t1), typeRegistry.getParticleTypeID(t2));
}

void PotentialRegistry::configure() {
    namespace coll = readdy::util::collections;
    using pair = util::particle_type_pair;
    using pot1 = potentials::PotentialOrder1;
    using pot1_ptr = std::unique_ptr<potentials::PotentialOrder1>;
    using pot2_ptr = std::unique_ptr<potentials::PotentialOrder2>;
    using pot2 = potentials::PotentialOrder2;
    potentialO1Registry.clear();
    potentialO2Registry.clear();

    coll::for_each_value(potentialO1RegistryInternal, [&](const particle_type_type type, const pot1_ptr &ptr) {
        ptr->configureForType(&typeRegistry, type);
        (potentialO1Registry)[type].push_back(ptr.get());
    });
    coll::for_each_value(potentialO2RegistryInternal, [&](const pair &type, const pot2_ptr &ptr) {
        ptr->configureForTypes(&typeRegistry, std::get<0>(type), std::get<1>(type));
        (potentialO2Registry)[type].push_back(ptr.get());
    });
    coll::for_each_value(potentialO1RegistryExternal, [&](const particle_type_type type, pot1 *ptr) {
        ptr->configureForType(&typeRegistry, type);
        (potentialO1Registry)[type].push_back(ptr);
    });
    coll::for_each_value(potentialO2RegistryExternal, [&](const pair &type, pot2 *ptr) {
        ptr->configureForTypes(&typeRegistry, std::get<0>(type), std::get<1>(type));
        (potentialO2Registry)[type].push_back(ptr);
    });
}

void PotentialRegistry::debug_output() const {
    auto find_pot_name = [this](particle_type_type type) -> const std::string {
        for (auto &&t : typeRegistry.getTypeMapping()) {
            if (t.second == type) return t.first;
        }
        return "";
    };
    if (!getAllOrder1Potentials().empty()) {
        log::debug(" - potentials of order 1:");
        for (const auto &types : getAllOrder1Potentials()) {
            log::debug("     * for type {}", find_pot_name(types.first));
            for (auto pot : types.second) {
                log::debug("         * {}", pot->describe());
            }
        }
    }
    if (!getAllOrder2Potentials().empty()) {
        log::debug(" - potentials of order 2:");
        for (const auto &types : getAllOrder2Potentials()) {
            log::debug("     * for types {} and {}", find_pot_name(std::get<0>(types.first)),
                       find_pot_name(std::get<1>(types.first)));
            for (auto pot : types.second) {
                log::debug("         * {}", pot->describe());
            }
        }
    }
}
}
}
}

