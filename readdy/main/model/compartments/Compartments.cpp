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
 * @file Compartments.cpp
 * @brief << brief description >>
 * @author chrisfroe
 * @date 13.01.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <readdy/model/compartments/Compartments.h>
#include <readdy/common/logging.h>

namespace readdy {
namespace model {
namespace compartments {

short Compartment::counter = 0;

Sphere::Sphere(const Compartment::conversion_map &conversions, const std::string &uniqueName, const Vec3 &origin,
               const scalar radius, const bool largerOrLess)
        : Compartment(conversions, getCompartmentTypeName<Sphere>(), uniqueName), radius(radius), radiusSquared(radius * radius),
          largerOrLess(largerOrLess), origin(origin) {}

Plane::Plane(const Compartment::conversion_map &conversions, const std::string &uniqueName, const Vec3 &normalCoefficients,
             const scalar distance, const bool largerOrLess)
        : Compartment(conversions, getCompartmentTypeName<Plane>(), uniqueName), normalCoefficients(normalCoefficients), distanceFromOrigin(distance),
          largerOrLess(largerOrLess) {
    const auto normSquared = normalCoefficients * normalCoefficients;
    if (std::abs(normSquared - 1) > 0.0001) {
        throw std::invalid_argument("Plane coefficients not sufficiently normalized. Make sure that coefficients and "
                                            "distance are according to the Hesse normal form");
    }
}

}
}
}