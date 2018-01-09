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

    virtual scalar getRelevantLengthScale() const noexcept override {
        return std::min(extent[0], std::min(extent[1], extent[2]));
    }

    virtual scalar getMaximalForce(scalar kbt) const noexcept override {
        return 0;
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

    void calculateForceAndEnergy(Vec3 &force, scalar &energy, const Vec3 &position) const override {
        energy += calculateEnergy(position);
        calculateForce(force, position);
    }

    std::string describe() const override;

    std::string type() const override;

protected:
    const Vec3 origin, extent, min, max;
    const scalar forceConstant;
};

// @todo modify this, so that you can choose whether the sphere keeps particles in or out
class SphereIn : public PotentialOrder1 {
    using super = PotentialOrder1;
public:
    SphereIn(particle_type_type particleType, scalar forceConstant, const Vec3& origin, scalar radius)
            : super(particleType), origin(origin), radius(radius), forceConstant(forceConstant) {}

    virtual scalar getRelevantLengthScale() const noexcept override {
        return radius;
    }

    virtual scalar getMaximalForce(scalar kbt) const noexcept override {
        return 0;
    }

    scalar calculateEnergy(const Vec3 &position) const override {
        auto difference = position - origin;
        scalar distanceFromOrigin = difference.norm();
        scalar distanceFromSphere = distanceFromOrigin - radius;
        auto energy = static_cast<scalar>(0.);
        if (distanceFromSphere > 0) {
            energy = static_cast<scalar>(0.5 * forceConstant * distanceFromSphere * distanceFromSphere);
        }
        return energy;
    }

    void calculateForce(Vec3 &force, const Vec3 &position) const override {
        auto difference = position - origin;
        scalar distanceFromOrigin = difference.norm();
        scalar distanceFromSphere = distanceFromOrigin - radius;
        if (distanceFromSphere > 0) {
            force += -1 * forceConstant * distanceFromSphere * difference / distanceFromOrigin;
        }
    }

    void calculateForceAndEnergy(Vec3 &force, scalar &energy, const Vec3 &position) const override {
        auto difference = position - origin;
        scalar distanceFromOrigin = difference.norm();
        scalar distanceFromSphere = distanceFromOrigin - radius;
        if (distanceFromSphere > 0) {
            energy += 0.5 * forceConstant * distanceFromSphere * distanceFromSphere;
            force += -1 * forceConstant * distanceFromSphere * difference / distanceFromOrigin;
        }
    }

    std::string describe() const override;

    std::string type() const override;

protected:
    const Vec3 origin;
    const scalar radius, forceConstant;
};

class SphereOut : public PotentialOrder1 {
    using super = PotentialOrder1;
public:
    SphereOut(particle_type_type particleType, scalar forceConstant, const Vec3& origin, scalar radius)
            : super(particleType), forceConstant(forceConstant), origin(origin), radius(radius) {}

    virtual scalar getRelevantLengthScale() const noexcept override {
        return radius;
    }

    virtual scalar getMaximalForce(scalar kbt) const noexcept override{
        return 0;
    }

    scalar calculateEnergy(const Vec3 &position) const override {
        auto difference = position - origin;
        scalar distanceFromOrigin = difference.norm();
        scalar distanceFromSphere = distanceFromOrigin - radius;
        auto energy = static_cast<scalar>(0.);
        if (distanceFromSphere < 0) {
            energy = static_cast<scalar>(0.5) * forceConstant * distanceFromSphere * distanceFromSphere;
        }
        return energy;
    }

    void calculateForce(Vec3 &force, const Vec3 &position) const override {
        auto difference = position - origin;
        scalar distanceFromOrigin = difference.norm();
        scalar distanceFromSphere = distanceFromOrigin - radius;
        if (distanceFromSphere < 0) {
            force += -1 * forceConstant * distanceFromSphere * difference / distanceFromOrigin;
        }
    }

    void calculateForceAndEnergy(Vec3 &force, scalar &energy, const Vec3 &position) const override {
        auto difference = position - origin;
        scalar distanceFromOrigin = difference.norm();
        scalar distanceFromSphere = distanceFromOrigin - radius;
        if (distanceFromSphere < 0) {
            energy += 0.5 * forceConstant * distanceFromSphere * distanceFromSphere;
            force += -1 * forceConstant * distanceFromSphere * difference / distanceFromOrigin;
        }
    }

    std::string describe() const override;

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

    virtual readdy::scalar getRelevantLengthScale() const noexcept override {
        return width;
    }

    virtual readdy::scalar getMaximalForce(readdy::scalar kbt) const noexcept override {
        return static_cast<scalar>(2. * height / width);
    }

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

    void calculateForceAndEnergy(Vec3 &force, readdy::scalar &energy, const Vec3 &position) const override {
        const auto difference = position - origin;
        const auto distance = difference.norm();
        if (distance < r1) {
            // nothing happens
        } else if (r4 <= distance) {
            // nothing happens
        } else if (r1 <= distance && distance < r2) {
            force += - effectiveForceConstant * (distance - radius + width) * difference / distance;
            energy += 0.5 * effectiveForceConstant * std::pow(distance - radius + width, 2.);
        } else if (r2 <= distance && distance < r3) {
            force += effectiveForceConstant * (distance - radius) * difference / distance;
            energy += height - 0.5 * effectiveForceConstant * std::pow(distance - radius, 2.);
        } else if (r3 <= distance && distance < r4) {
            force += - effectiveForceConstant * (distance - radius - width) * difference / distance;
            energy += 0.5 * effectiveForceConstant * std::pow(distance - radius - width, 2.);
        } else {
            throw std::runtime_error("Critical error in calculateForceAndEnergy of spherical barrier");
        }
    }

    std::string describe() const override;

    std::string type() const override;

protected:
    const Vec3 origin;
    const readdy::scalar radius, height, width, r1, r2, r3, r4, effectiveForceConstant;
};

template<typename T>
const std::string getPotentialName(typename std::enable_if<std::is_base_of<Box, T>::value>::type * = 0) {
    return "Box";
}
template<typename T>
const std::string getPotentialName(typename std::enable_if<std::is_base_of<SphereIn, T>::value>::type* = 0) {
    return "SphereIn";
}
template <typename T>
const std::string getPotentialName(typename std::enable_if<std::is_base_of<SphereOut, T>::value>::type* = 0) {
    return "SphereOut";
}
template <typename T>
const std::string getPotentialName(typename std::enable_if<std::is_base_of<SphericalBarrier, T>::value>::type * = 0) {
    return "SphericalBarrier";
}

NAMESPACE_END(potentials)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
