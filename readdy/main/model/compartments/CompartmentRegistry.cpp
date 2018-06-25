/********************************************************************
 * Copyright © 2017 Computational Molecular Biology Group,          * 
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
 * @file CompartmentRegistry.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 20.09.17
 * @copyright GNU Lesser General Public License v3.0
 */

#include <readdy/model/compartments/CompartmentRegistry.h>
#include <readdy/model/compartments/Compartments.h>
#include <readdy/model/_internal/Util.h>

namespace readdy {
namespace model {
namespace compartments {


Compartment::id_type
CompartmentRegistry::addSphere(const Compartment::conversion_map &conversions, const std::string &uniqueName,
                               const Vec3 &origin, scalar radius, bool largerOrLess) {
    _compartments.emplace_back(std::make_shared<Sphere>(conversions, uniqueName, origin, radius, largerOrLess));
    return _compartments.back()->getId();
}

Compartment::id_type
CompartmentRegistry::addPlane(const Compartment::conversion_map &conversions, const std::string &uniqueName,
                              const Vec3 &normalCoefficients, scalar distance, bool largerOrLess) {
    _compartments.emplace_back(std::make_shared<Plane>(conversions, uniqueName, normalCoefficients, distance,
                                                       largerOrLess));
    return _compartments.back()->getId();
}

CompartmentRegistry::CompartmentRegistry(const ParticleTypeRegistry &types) : _types(types) {}

}
}
}
