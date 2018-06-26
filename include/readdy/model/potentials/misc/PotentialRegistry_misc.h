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

inline void PotentialRegistry::remove(const PotentialId handle) {
    for (auto &entry : _ownPotentialsO1) {
        entry.second.erase(std::remove_if(entry.second.begin(), entry.second.end(),
                                          [=](const std::shared_ptr<potentials::PotentialOrder1> &p) -> bool {
                                              return handle == p->getId();
                                          }
        ), entry.second.end());
    }
    for (auto &entry : _ownPotentialsP2) {
        entry.second.erase(std::remove_if(entry.second.begin(), entry.second.end(),
                                          [=](const std::shared_ptr<potentials::PotentialOrder2> &p) -> bool {
                                              return handle == p->getId();
                                          }
        ), entry.second.end());
    }
    for (auto &entry : _potentialsO1) {
        entry.second.erase(std::remove_if(entry.second.begin(), entry.second.end(),
                                          [=](potentials::PotentialOrder1 *p) -> bool {
                                              return handle == p->getId();
                                          }
        ), entry.second.end());
    }
    for (auto &entry : _potentialsO2) {
        entry.second.erase(std::remove_if(entry.second.begin(), entry.second.end(),
                                          [=](potentials::PotentialOrder2 *p) -> bool {
                                              return handle == p->getId();
                                          }
        ), entry.second.end());
    }
}

inline std::string PotentialRegistry::describe() const {
    namespace rus = readdy::util::str;
    auto find_pot_name = [this](ParticleTypeId type) -> const std::string {
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
            description += fmt::format(R"(     * for types "{}" and "{}"{})", find_pot_name(std::get<0>(types.first)),
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
