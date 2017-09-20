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
 *
 * @file CompartmentFactory.h
 * @brief Create-methods and dispatchers for compartments
 * @author chrisfroe
 * @date 18.01.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include "Compartments.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(compartments)

class CompartmentFactory {
public:
    template<typename T, typename... Args>
    std::shared_ptr<T> createCompartment(const Compartment::conversion_map &convMap, Args &&... args) const {
        return std::shared_ptr<T>(get_dispatcher<T, Args...>::impl(this, convMap, std::forward<Args>(args)...));
    }

protected:
    virtual Plane *
    createPlane(const Compartment::conversion_map &convMap, const std::string &uniqueName, const Vec3 &coefficients,
                const scalar distance, const bool largerOrLess) const {
        return new Plane(convMap, uniqueName, coefficients, distance, largerOrLess);
    }

    virtual Sphere *createSphere(const Compartment::conversion_map &convMap, const std::string &uniqueName, const Vec3 &origin,
                                 const scalar radius, const bool largerOrLess) const {
        return new Sphere(convMap, uniqueName, origin, radius, largerOrLess);
    }

    template<typename T, typename... Args>
    struct get_dispatcher;

    // default dispatcher only invokes normal constructor
    template<typename T, typename... Args>
    struct get_dispatcher {
        static T *impl(const CompartmentFactory *self, const Compartment::conversion_map &convMap, Args &&... args) {
            return new T(convMap, std::forward<Args>(args)...);
        };
    };

};

// specialized dispatcher methods will call the kernel-specific 'create...' methods
READDY_CREATE_COMPARTMENT_FACTORY_DISPATCHER(Plane)

READDY_CREATE_COMPARTMENT_FACTORY_DISPATCHER(Sphere)

NAMESPACE_END(compartments)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
