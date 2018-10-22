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
 * This header file contains the declarations of all possibly available order 1 potentials. Currently:
 *   - Box potential
 *   - Sphere potential
 *
 * @file PotentialsOrder1.h
 * @brief This header file contains the declarations of all possibly available order 1 potentials.
 * @author clonker
 * @author chrisfroe
 * @date 15.06.16
 */

#pragma once

#include <ostream>
#include "PotentialOrder1.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(potentials)

// @todo allow box exclusion
class Box : public PotentialOrder1 {
    using super = PotentialOrder1;
public:
    Box(particle_type_type particleType, scalar forceConstant, const Vec3& origin, const Vec3& extent);

    const Vec3 &getOrigin() const {
        return origin;
    }

    const Vec3 &getExtent() const {
        return extent;
    }

    scalar getForceConstant() const {
        return forceConstant;
    }

    scalar calculateEnergy(const Vec3 &position) const override {
        scalar energy = 0;

        for (auto i = 0; i < 3; ++i) {
            if (position[i] < min[i] || position[i] > max[i]) {
                if (position[i] < min[i]) {
                    energy += 0.5 * forceConstant * (position[i] - min[i]) * (position[i] - min[i]);
                } else {
                    energy += 0.5 * forceConstant * (position[i] - max[i]) * (position[i] - max[i]);
                }
            }
        }

        return energy;
    }

    void calculateForce(Vec3 &force, const Vec3 &position) const override {
        for (auto i = 0; i < 3; i++) {
            if (position[i] < min[i] || position[i] > max[i]) {
                if (position[i] < min[i]) {
                    force[i] += -1 * forceConstant * (position[i] - min[i]);
                } else {
                    force[i] += -1 * forceConstant * (position[i] - max[i]);
                }
            }
        }
    }

    std::string describe() const override;

    std::string type() const override;

protected:
    const Vec3 origin, extent, min, max;
    const scalar forceConstant;
};

template<bool inclusion>
class Sphere : public PotentialOrder1 {
    using super = PotentialOrder1;
public:
    Sphere(particle_type_type particleType, scalar forceConstant, const Vec3& origin, scalar radius)
            : super(particleType), origin(origin), radius(radius), forceConstant(forceConstant) {}

    scalar calculateEnergy(const Vec3 &position) const override {
        auto difference = position - origin;
        scalar distanceFromOrigin = difference.norm();
        scalar distanceFromSphere = distanceFromOrigin - radius;
        if constexpr (inclusion) {
            if (distanceFromSphere > 0) {
                return static_cast<scalar>(0.5) * forceConstant * distanceFromSphere * distanceFromSphere;
            }
        } else if constexpr (!inclusion) {
            if (distanceFromSphere < 0) {
                return static_cast<scalar>(0.5) * forceConstant * distanceFromSphere * distanceFromSphere;
            }
        }
        return static_cast<scalar>(0.);
    }

    void calculateForce(Vec3 &force, const Vec3 &position) const override {
        auto difference = position - origin;
        scalar distanceFromOrigin = difference.norm();
        scalar distanceFromSphere = distanceFromOrigin - radius;
        if constexpr (inclusion) {
            if (distanceFromSphere > 0) {
                force += -1 * forceConstant * distanceFromSphere * difference / distanceFromOrigin;
            }
        } else if constexpr (!inclusion) {
            if (distanceFromSphere < 0) {
                force += -1 * forceConstant * distanceFromSphere * difference / distanceFromOrigin;
            }
        }
    }

    std::string describe() const override {
        std::string inOrOut = inclusion ? "inclusion" : "exclusion";
        return fmt::format("Spherical {} potential with origin={}, radius={}, and force constant k={}",
                           inOrOut, origin, radius, forceConstant);
    }

    std::string type() const override;

protected:
    const Vec3 origin;
    const scalar radius, forceConstant;
};

/**
 * A potential that forms a concentric barrier at a certain radius around a given origin. It is given a height (in terms of energy)
 * and a width. Note that the height can also be negative, then this potential acts as a 'sticky' sphere. The potential consists
 * of harmonic snippets, such that the energy landscape is continuous and differentiable, the force is only continuous and not differentiable.
 */
class SphericalBarrier : public PotentialOrder1 {
    using super = PotentialOrder1;
public:
    SphericalBarrier(particle_type_type particleType, scalar height, scalar width, const Vec3 &origin, scalar radius);

    readdy::scalar calculateEnergy(const Vec3 &position) const override {
        const auto difference = position - origin;
        const auto distance = difference.norm();
        if (distance < r1) {
            return static_cast<scalar>(0.);
        }
        if (r4 <= distance) {
            return static_cast<scalar>(0.);
        }
        if (r1 <= distance && distance < r2) {
            return static_cast<scalar>(0.5) * effectiveForceConstant * std::pow(distance - radius + width, static_cast<scalar>(2.));
        }
        if (r2 <= distance && distance < r3) {
            return height - static_cast<scalar>(0.5) * effectiveForceConstant * std::pow(distance - radius, static_cast<scalar>(2.));
        }
        if (r3 <= distance && distance < r4) {
            return static_cast<scalar>(0.5) * effectiveForceConstant * std::pow(distance - radius - width, static_cast<scalar>(2.));
        }
        throw std::runtime_error("Critical error in calculateEnergy of spherical barrier");
    }

    void calculateForce(Vec3 &force, const Vec3 &position) const override {
        const auto difference = position - origin;
        const auto distance = difference.norm();
        if (distance < r1) {
            // nothing happens
        } else if (r4 <= distance) {
            // nothing happens
        } else if (r1 <= distance && distance < r2) {
            force += - effectiveForceConstant * (distance - radius + width) * difference / distance;
        } else if (r2 <= distance && distance < r3) {
            force += effectiveForceConstant * (distance - radius) * difference / distance;
        } else if (r3 <= distance && distance < r4) {
            force += - effectiveForceConstant * (distance - radius - width) * difference / distance;
        } else {
            throw std::runtime_error("Critical error in calculateForce of spherical barrier");
        }
    }

    std::string describe() const override;

    std::string type() const override;

protected:
    const Vec3 origin;
    const readdy::scalar radius, height, width, r1, r2, r3, r4, effectiveForceConstant;
};

template<bool inclusion>
class Cylinder : public PotentialOrder1 {
    using super = PotentialOrder1;
public:
    Cylinder(particle_type_type particleType, scalar forceConstant, const Vec3 &origin, const Vec3 &normal,
             scalar radius) : super(particleType), forceConstant(forceConstant),
                              origin(origin), normal(normal / normal.norm()), radius(radius) {};

    scalar calculateEnergy(const Vec3 &position) const override {
        Vec3 pos = position - origin;
        Vec3 perpendicular = pos - (pos * normal) * normal;
        scalar distanceFromAxis = perpendicular.norm();
        if constexpr (inclusion) {
            if (distanceFromAxis > radius) {
                return static_cast<scalar>(0.5) * forceConstant *
                       std::pow(distanceFromAxis - radius, static_cast<scalar>(2));
            }
        } else if constexpr (!inclusion) {
            if (distanceFromAxis < radius) {
                return static_cast<scalar>(0.5) * forceConstant *
                       std::pow(distanceFromAxis - radius, static_cast<scalar>(2));
            }
        }
        return static_cast<scalar>(0.);
    }

    void calculateForce(Vec3 &force, const Vec3 &position) const override {
        Vec3 pos = position - origin;
        Vec3 perpendicular = pos - (pos * normal) * normal;
        scalar distanceFromAxis = perpendicular.norm();
        if constexpr (inclusion) {
            if (distanceFromAxis > radius) {
                force += -forceConstant * (distanceFromAxis - radius) * perpendicular / distanceFromAxis;
            }
        } else if constexpr (!inclusion) {
            if (distanceFromAxis < radius) {
                force += -forceConstant * (distanceFromAxis - radius) * perpendicular / distanceFromAxis;
            }
        }
    }

    std::string describe() const override {
        std::string inOrOut = inclusion ? "inclusion" : "exclusion";
        return fmt::format("Cylindrical {} potential with force constant={}, origin={}, normal={}, and radius={}",
                           inOrOut, forceConstant, origin, normal, radius);
    }

    std::string type() const override;

protected:
    const Vec3 origin;
    const Vec3 normal;
    const scalar forceConstant, radius;
};

template<typename T>
const std::string getPotentialName(typename std::enable_if<std::is_base_of<Box, T>::value>::type * = 0) {
    return "Box";
}
template<typename T>
const std::string getPotentialName(typename std::enable_if<std::is_base_of<Sphere<true>, T>::value>::type* = 0) {
    return "SphereInclusion";
}
template<typename T>
const std::string getPotentialName(typename std::enable_if<std::is_base_of<Sphere<false>, T>::value>::type* = 0) {
    return "SphereExclusion";
}
template <typename T>
const std::string getPotentialName(typename std::enable_if<std::is_base_of<SphericalBarrier, T>::value>::type * = 0) {
    return "SphericalBarrier";
}
template <typename T>
const std::string getPotentialName(typename std::enable_if<std::is_base_of<Cylinder<true>, T>::value>::type * = 0) {
    return "CylinderInclusion";
}
template <typename T>
const std::string getPotentialName(typename std::enable_if<std::is_base_of<Cylinder<false>, T>::value>::type * = 0) {
    return "CylinderExclusion";
}


NAMESPACE_END(potentials)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
