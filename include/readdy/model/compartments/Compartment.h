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
 * @file Compartment.h
 * @brief << brief description >>
 * @author chrisfroe
 * @date 12.01.17
 * @copyright GNU Lesser General Public License v3.0
 */

#ifndef READDY_MAIN_COMPARTMENT_H
#define READDY_MAIN_COMPARTMENT_H

#include <readdy/model/Vec3.h>
#include <readdy/model/Particle.h>
#include <unordered_map>

namespace readdy {
namespace model {
namespace compartments {

class Compartment {
public:
    using id_t = short;
    using particleType_t = readdy::model::Particle::type_type;

    Compartment(const std::unordered_map<particleType_t, particleType_t> &conversions, const std::string &typeName, const std::string &uniqueName)
            : typeName(typeName), uniqueName(uniqueName), conversions(conversions), id(counter++) {}

    virtual const bool isContained(const Vec3 &position) const = 0;

    const std::unordered_map<particleType_t, particleType_t> &getConversions() const {
        return conversions;
    }

    const std::string &getTypeName() const {
        return typeName;
    }

    const std::string &getUniqueName() const {
        return uniqueName;
    }

    const short getId() const {
        return id;
    }

private:
    static id_t counter;

    const std::string typeName;
    const std::string uniqueName;
    const id_t id;
    const std::unordered_map<particleType_t, particleType_t> conversions;
};

}
}
}

#endif //READDY_MAIN_COMPARTMENT_H
