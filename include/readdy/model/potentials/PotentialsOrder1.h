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

    const Vec3 &getOrigin() const;

    const Vec3 &getExtent() const;

    scalar getForceConstant() const;

    virtual scalar getRelevantLengthScale() const noexcept override;

    virtual scalar getMaximalForce(scalar kbt) const noexcept override;

    scalar calculateEnergy(const Vec3 &position) const override;

    void calculateForce(Vec3 &force, const Vec3 &position) const override;

    void calculateForceAndEnergy(Vec3 &force, scalar &energy, const Vec3 &position) const override;

    std::string describe() const override;

protected:
    const Vec3 origin, extent, min, max;
    const scalar forceConstant;
};

// @todo modify this, so that you can choose whether the sphere keeps particles in or out
class SphereIn : public PotentialOrder1 {
    using super = PotentialOrder1;
public:
    SphereIn(particle_type_type particleType, scalar forceConstant, const Vec3& origin, scalar radius);

    virtual scalar getRelevantLengthScale() const noexcept override;

    virtual scalar getMaximalForce(scalar kbt) const noexcept override;

    scalar calculateEnergy(const Vec3 &position) const override;

    void calculateForce(Vec3 &force, const Vec3 &position) const override;

    void calculateForceAndEnergy(Vec3 &force, scalar &energy, const Vec3 &position) const override;

    std::string describe() const override;

protected:
    const Vec3 origin;
    const scalar radius, forceConstant;
};

class SphereOut : public PotentialOrder1 {
    using super = PotentialOrder1;
public:
    SphereOut(particle_type_type particleType, scalar forceConstant, const Vec3& origin, scalar radius);

    virtual scalar getRelevantLengthScale() const noexcept override;

    virtual scalar getMaximalForce(scalar kbt) const noexcept override;

    scalar calculateEnergy(const Vec3 &position) const override;

    void calculateForce(Vec3 &force, const Vec3 &position) const override;

    void calculateForceAndEnergy(Vec3 &force, scalar &energy, const Vec3 &position) const override;

    std::string describe() const override;
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

    virtual readdy::scalar getRelevantLengthScale() const noexcept override;

    virtual readdy::scalar getMaximalForce(readdy::scalar kbt) const noexcept override;

    readdy::scalar calculateEnergy(const Vec3 &position) const override;

    void calculateForce(Vec3 &force, const Vec3 &position) const override;

    void calculateForceAndEnergy(Vec3 &force, readdy::scalar &energy, const Vec3 &position) const override;

    std::string describe() const override;

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
