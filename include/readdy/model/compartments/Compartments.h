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

#include "Compartment.h"

namespace readdy::model::compartments {

class Sphere : public Compartment {
public:
    Sphere(const conversion_map &conversions, const std::string &uniqueName, const Vec3 &origin,
           scalar radius, bool largerOrLess);

    bool isContained(const Vec3 &position) const override {
        const auto delta = position - origin;
        const auto distanceSquared = delta * delta;
        if (largerOrLess) {
            return distanceSquared > radiusSquared;
        }
        return distanceSquared < radiusSquared;
    }

protected:
    const Vec3 origin;
    const scalar radius;
    const scalar radiusSquared;
    const bool largerOrLess;
};


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

class Capsule : public Compartment {
public:
    Capsule(const conversion_map &conversions, const std::string &uniqueName, Vec3 center, Vec3 direction,
            scalar length, scalar radius, bool inside)
        : Compartment(conversions, "Capsule", uniqueName), _center(center),
          _direction(direction / direction.norm()), _radius(radius), _length(length), _inside(inside) {

    }

    template<typename T>
    auto differentSign(T t1, T t2) const {
        return (t1 >= 0 && t2 < 0) || (t1 < 0 && t2 >= 0);
    }

    auto closestCircleCenter(const Vec3 &position) const {
        // we define the line x(l) = x_0 + l*v, where v is the normalized direction vector
        // then for given position the lambda is (x_0 - p)*v*v:
        auto v = (((position - _center) * _direction) * _direction);
        auto lambda = v.norm();
        // check if any of the signs in v differ from direction:
        if (differentSign(v.x, _direction.x) || differentSign(v.y, _direction.y) || differentSign(v.z, _direction.z)) {
            lambda = -lambda;
        }
        // clamp lambda to half length of capsule
        lambda = std::clamp(lambda, -_length / 2, _length / 2);
        // now we obtain point corresponding to lambda
        auto circleCenter = _center + lambda * _direction;
        return circleCenter;
    }

    auto distToCapsuleSquared(const Vec3 &position) const {
        auto cc = closestCircleCenter(position);
        // and check if point is inside circle
        auto distToCCSquared = (cc - position).normSquared();
        return distToCCSquared;
    }

    bool isContained(const Vec3 &position) const override {
        auto insideSphere = distToCapsuleSquared(position) <= _radius*_radius;
        return _inside == insideSphere;
    }

    const auto &center() const { return _center; }
    const auto &direction() const { return _direction; }
    auto radius() const { return _radius; }
    auto length() const { return _length; }
    auto inside() const { return _inside; }

protected:
    Vec3 _center, _direction;
    scalar _radius, _length;
    bool _inside;
};

template<typename T>
std::string getCompartmentTypeName(typename std::enable_if<std::is_base_of<Sphere, T>::value>::type * = 0) {
    return "Sphere";
}

template<typename T>
std::string getCompartmentTypeName(typename std::enable_if<std::is_base_of<Plane, T>::value>::type * = 0) {
    return "Plane";
}

}
