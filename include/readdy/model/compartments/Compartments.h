/********************************************************************
 * Copyright © 2018 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * Redistribution and use in source and binary forms, with or       *
 * without modification, are permitted provided that the            *
 * following conditions are met:                                    *
 *  1. Redistributions of source code must retain the above         *
 *     copyright notice, this list of conditions and the            *
 *     following disclaimer.                                        *
 *  2. Redistributions in binary form must reproduce the above      *
 *     copyright notice, this list of conditions and the following  *
 *     disclaimer in the documentation and/or other materials       *
 *     provided with the distribution.                              *
 *  3. Neither the name of the copyright holder nor the names of    *
 *     its contributors may be used to endorse or promote products  *
 *     derived from this software without specific                  *
 *     prior written permission.                                    *
 *                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND           *
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,      *
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF         *
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE         *
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR            *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,     *
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,         *
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER *
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,      *
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)    *
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF      *
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
 ********************************************************************/


/**
 * << detailed description >>
 *
 * @file Compartments.h
 * @brief << brief description >>
 * @author chrisfroe
 * @date 13.01.17
 * @copyright BSD-3
 */

#pragma once

#include "../geometry.h"
#include "Compartment.h"

namespace readdy::model::compartments {


template<typename Geometry>
class GeometryCompartment : public Compartment {
public:
    GeometryCompartment(const conversion_map &conversions, const std::string &uniqueName,
                        const Geometry &geometry, bool inside)
    : Compartment(conversions, std::string(Geometry::name), uniqueName), geometry(geometry), inside(inside) {}

    [[nodiscard]] bool isContained(const Vec3 &position) const override {
        if(inside) {
            return geometry.template contains<true>(position);
        } else {
            return geometry.template contains<false>(position);
        }
    }

private:
    Geometry geometry;
    bool inside;
};

using Sphere = GeometryCompartment<geometry::Sphere<scalar>>;
using Box = GeometryCompartment<geometry::Box<scalar>>;
using Capsule = GeometryCompartment<geometry::Capsule<scalar>>;

class Plane : public Compartment {
public:
    Plane(const conversion_map &conversions, const std::string &uniqueName, const Vec3 &normalCoefficients,
          scalar distance, bool largerOrLess);

    bool isContained(const Vec3 &position) const override {
        const scalar distanceFromPlane = position * normalCoefficients - distanceFromOrigin;
        if (largerOrLess) {
            return distanceFromPlane > 0;
        }
        return distanceFromPlane < 0;
    }

protected:
    const Vec3 normalCoefficients;
    const scalar distanceFromOrigin;
    const bool largerOrLess;
};

}
