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
#pragma once

#include <unordered_map>
#include <utility>

#include <readdy/model/Particle.h>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(compartments)

class Compartment {
public:
    using id_type = short;
    using label_conversion_map = std::unordered_map<std::string, std::string>;
    using conversion_map = std::unordered_map<particle_type_type, particle_type_type>;

    Compartment(conversion_map conversions, std::string typeName,
                std::string uniqueName)
            : typeName(std::move(typeName)), uniqueName(std::move(uniqueName)),
              conversions(std::move(conversions)), _id(counter++) {}

    virtual const bool isContained(const Vec3 &position) const = 0;

    const conversion_map &getConversions() const {
        return conversions;
    }

    const std::string &getTypeName() const {
        return typeName;
    }

    const std::string &getUniqueName() const {
        return uniqueName;
    }

    const id_type getId() const {
        return _id;
    }

protected:
    static id_type counter;

    std::string typeName;
    std::string uniqueName;
    id_type _id;
    conversion_map conversions;
};

NAMESPACE_END(compartments)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
