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

const Potential::id PotentialRegistry::add_external(potentials::PotentialOrder2 *potential) {
    const auto id = potential->getId();
    auto type1Id = typeRegistry.id_of(potential->particleType1);
    auto type2Id = typeRegistry.id_of(potential->particleType2);
    auto pp = std::tie(type1Id, type2Id);
    if (potentialO2RegistryExternal.find(pp) == potentialO2RegistryExternal.end()) {
        potentialO2RegistryExternal.emplace(pp, pot_ptr_vec2_external());
    }
    potentialO2RegistryExternal[pp].push_back(potential);
    return id;
}

const Potential::id PotentialRegistry::add_external(potentials::PotentialOrder1 *potential) {
    auto typeId = typeRegistry.id_of(potential->particleType);
    if (potentialO1RegistryExternal.find(typeId) == potentialO1RegistryExternal.end()) {
        potentialO1RegistryExternal.emplace(std::make_pair(typeId, pot_ptr_vec1_external()));
    }
    potentialO1RegistryExternal[typeId].push_back(potential);
    return potential->getId();
}

void PotentialRegistry::remove(const Potential::id handle) {
    for (auto &entry : potentialO1RegistryInternal) {
        entry.second.erase(std::remove_if(entry.second.begin(), entry.second.end(),
                                        [=](const std::unique_ptr<potentials::PotentialOrder1> &p) -> bool {
                                            return handle == p->getId();
                                        }
        ), entry.second.end());
    }
    for (auto &entry : potentialO2RegistryInternal) {
        entry.second.erase(std::remove_if(entry.second.begin(), entry.second.end(),
                                        [=](const std::unique_ptr<potentials::PotentialOrder2> &p) -> bool {
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

const PotentialRegistry::potentials_o1& PotentialRegistry::potentials_of(const particle_type_type type) const {
    auto it = potentialO1Registry.find(type);
    if(it != potentialO1Registry.end()) {
        return it->second;
    }
    return defaultPotentialsO1;
}

const PotentialRegistry::potential_o1_registry& PotentialRegistry::potentials_order1() const {
    return potentialO1Registry;
}

const std::vector<PotentialOrder2 *> &
PotentialRegistry::potentials_of(const particle_type_type t1, const particle_type_type t2) const {
    auto it = potentialO2Registry.find(std::tie(t1, t2));
    if(it != potentialO2Registry.end()) {
        return it->second;
    }
    return defaultPotentialsO2;
}

const PotentialRegistry::potential_o2_registry &PotentialRegistry::potentials_order2() const {
    return potentialO2Registry;
}

const PotentialRegistry::potentials_o1& PotentialRegistry::potentials_of(const std::string &type) const {
    return potentials_of(typeRegistry.id_of(type));
}

const PotentialRegistry::potentials_o2 &
PotentialRegistry::potentials_of(const std::string &t1, const std::string &t2) const {
    return potentials_of(typeRegistry.id_of(t1), typeRegistry.id_of(t2));
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
        for (auto &&t : typeRegistry.type_mapping()) {
            if (t.second == type) return t.first;
        }
        return "";
    };
    if (!potentials_order1().empty()) {
        log::debug(" - potentials of order 1:");
        for (const auto &types : potentials_order1()) {
            log::debug("     * for type {}", find_pot_name(types.first));
            for (auto pot : types.second) {
                log::debug("         * {}", pot->describe());
            }
        }
    }
    if (!potentials_order2().empty()) {
        log::debug(" - potentials of order 2:");
        for (const auto &types : potentials_order2()) {
            log::debug("     * for types {} and {}", find_pot_name(std::get<0>(types.first)),
                       find_pot_name(std::get<1>(types.first)));
            for (auto pot : types.second) {
                log::debug("         * {}", pot->describe());
            }
        }
    }
}

const Potential::id PotentialRegistry::add(std::unique_ptr<PotentialOrder1> potential) {
    const auto id = potential->getId();
    auto typeId = typeRegistry.id_of(potential->particleType);
    potentialO1RegistryInternal[typeId].push_back(std::move(potential));
    return id;
}

const Potential::id PotentialRegistry::add(std::unique_ptr<PotentialOrder2> potential) {
    const auto id = potential->getId();
    auto type1Id = typeRegistry.id_of(potential->particleType1);
    auto type2Id = typeRegistry.id_of(potential->particleType2);
    potentialO2RegistryInternal[std::tie(type1Id, type2Id)].push_back(std::move(potential));
    return id;
}
}
}
}

