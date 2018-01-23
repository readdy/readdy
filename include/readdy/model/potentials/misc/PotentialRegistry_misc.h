/********************************************************************
 * Copyright © 2018 Computational Molecular Biology Group,          *
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
 * @file PotentialRegistry_misc.h
 * @brief << brief description >>
 * @author clonker
 * @date 1/11/18
 */


#pragma once

#include <readdy/common/Utils.h>
#include <readdy/common/string.h>
#include "../PotentialRegistry.h"

namespace readdy {
namespace model {
namespace potentials {

inline PotentialRegistry::id_type PotentialRegistry::addUserDefined(potentials::PotentialOrder2 *potential) {
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

inline PotentialRegistry::id_type PotentialRegistry::addUserDefined(potentials::PotentialOrder1 *potential) {
    auto typeId = potential->particleType();
    if (potentialO1RegistryExternal.find(typeId) == potentialO1RegistryExternal.end()) {
        potentialO1RegistryExternal.emplace(std::make_pair(typeId, pot_ptr_vec1_external()));
    }
    potentialO1RegistryExternal[typeId].push_back(potential);
    return potential->getId();
}

inline void PotentialRegistry::remove(const Potential::id_type handle) {
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

inline void PotentialRegistry::configure() {
    namespace coll = readdy::util::collections;
    using pair = util::particle_type_pair;
    using pot1 = potentials::PotentialOrder1;
    using pot1_ptr = std::shared_ptr<potentials::PotentialOrder1>;
    using pot2_ptr = std::shared_ptr<potentials::PotentialOrder2>;
    using pot2 = potentials::PotentialOrder2;
    potentialO1Registry.clear();
    potentialO2Registry.clear();
    _alternativeO2Registry.clear();

    coll::for_each_value(potentialO1RegistryInternal, [&](const particle_type_type type, const pot1_ptr &ptr) {
        (potentialO1Registry)[type].push_back(ptr.get());
    });
    coll::for_each_value(potentialO2RegistryInternal, [&](const pair &type, const pot2_ptr &ptr) {
        (potentialO2Registry)[type].push_back(ptr.get());
        _alternativeO2Registry[std::get<0>(type)][std::get<1>(type)].push_back(ptr.get());
        if(std::get<0>(type) != std::get<1>(type)) {
            _alternativeO2Registry[std::get<1>(type)][std::get<0>(type)].push_back(ptr.get());
        }
    });
    coll::for_each_value(potentialO1RegistryExternal, [&](const particle_type_type type, pot1 *ptr) {
        (potentialO1Registry)[type].push_back(ptr);
    });
    coll::for_each_value(potentialO2RegistryExternal, [&](const pair &type, pot2 *ptr) {
        (potentialO2Registry)[type].push_back(ptr);
        _alternativeO2Registry[std::get<0>(type)][std::get<1>(type)].push_back(ptr);
        if(std::get<0>(type) != std::get<1>(type)) {
            _alternativeO2Registry[std::get<1>(type)][std::get<0>(type)].push_back(ptr);
        }
    });
}

inline std::string PotentialRegistry::describe() const {
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
            description += fmt::format("     * for type \"{}\"{}", find_pot_name(types.first), rus::newline);
            for (auto pot : types.second) {
                description += fmt::format("         * {}{}", pot->describe(), rus::newline);
            }
        }
    }
    if (!potentialsOrder2().empty()) {
        description += fmt::format(" - potentials of order 2:{}", rus::newline);
        for (const auto &types : potentialsOrder2()) {
            description += fmt::format("     * for types \"{}\" and \"{}\"{}", find_pot_name(std::get<0>(types.first)),
                                       find_pot_name(std::get<1>(types.first)), rus::newline);
            for (auto pot : types.second) {
                description += fmt::format("         * {}{}", pot->describe(), rus::newline);
            }
        }
    }
    return description;
}

}
}
}
