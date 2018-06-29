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
 * @file CompartmentRegistry.h
 * @brief << brief description >>
 * @author clonker
 * @date 20.09.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <readdy/common/macros.h>
#include <readdy/model/ParticleTypeRegistry.h>
#include <readdy/model/_internal/Util.h>
#include "Compartment.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(compartments)

class CompartmentRegistry {
public:

    using CompartmentVector = std::vector<std::shared_ptr<readdy::model::compartments::Compartment>>;

    explicit CompartmentRegistry(const ParticleTypeRegistry &types);

    Compartment::id_type addSphere(const Compartment::conversion_map &conversions, const std::string &uniqueName,
                                   const Vec3 &origin, scalar radius, bool largerOrLess);

    Compartment::id_type addSphere(const Compartment::label_conversion_map &conversions, const std::string &uniqueName,
                                   const Vec3 &origin, scalar radius, bool largerOrLess) {
        return addSphere(_internal::util::transformTypesMap(conversions, _types.get()), uniqueName, origin, radius,
                         largerOrLess);
    }

    Compartment::id_type addPlane(const Compartment::conversion_map &conversions, const std::string &uniqueName,
                                  const Vec3 &normalCoefficients, scalar distance, bool largerOrLess);

    Compartment::id_type addPlane(const Compartment::label_conversion_map &conversions, const std::string &uniqueName,
                                  const Vec3 &normalCoefficients, scalar distance, bool largerOrLess) {
        return addPlane(_internal::util::transformTypesMap(conversions, _types.get()), uniqueName, normalCoefficients,
                        distance, largerOrLess);
    }


    const CompartmentVector &get() const {
        return _compartments;
    }

    CompartmentVector &get() {
        return _compartments;
    }

private:
     CompartmentVector _compartments;

    std::reference_wrapper<const ParticleTypeRegistry> _types;
};

NAMESPACE_END(compartments)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
