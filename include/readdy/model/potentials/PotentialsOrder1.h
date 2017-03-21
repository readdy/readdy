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
 *   - Cube potential
 *   - Sphere potential
 *
 * @file PotentialsOrder1.h
 * @brief This header file contains the declarations of all possibly available order 1 potentials.
 * @author clonker
 * @author chrisfroe
 * @date 15.06.16
 */

#include <ostream>
#include "PotentialOrder1.h"

#pragma once

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(potentials)

class Cube : public PotentialOrder1 {
    using super = PotentialOrder1;
public:
    Cube(const std::string& particleType, double forceConstant, const Vec3& origin, const Vec3& extent,
                  bool considerParticleRadius = true);

    const Vec3 &getOrigin() const;

    const Vec3 &getExtent() const;

    double getForceConstant() const;

    bool isConsiderParticleRadius() const;

    double getParticleRadius() const;

    virtual double getRelevantLengthScale() const noexcept override;

    virtual double getMaximalForce(double kbt) const noexcept override;

    double calculateEnergy(const Vec3 &position) const override;

    void calculateForce(Vec3 &force, const Vec3 &position) const override;

    void calculateForceAndEnergy(Vec3 &force, double &energy, const Vec3 &position) const override;

    void describe(std::ostream &os) const override;

protected:
    friend class readdy::model::KernelContext;

    void configureForType(const KernelContext *const ctx, const particle_type_type type) override;

    const Vec3 origin, extent, min, max;
    const double forceConstant;
    const bool considerParticleRadius;
    double particleRadius;
};

// @todo modify this, so that you can choose whether the sphere keeps particles in or out
class SphereIn : public PotentialOrder1 {
    using super = PotentialOrder1;
public:
    SphereIn(const std::string& particleType, double forceConstant, const Vec3& origin, double radius);

    virtual double getRelevantLengthScale() const noexcept override;

    virtual double getMaximalForce(double kbt) const noexcept override;

    double calculateEnergy(const Vec3 &position) const override;

    void calculateForce(Vec3 &force, const Vec3 &position) const override;

    void calculateForceAndEnergy(Vec3 &force, double &energy, const Vec3 &position) const override;

    void describe(std::ostream &os) const override;

protected:
    friend class readdy::model::KernelContext;

    void configureForType(const KernelContext *const ctx, const particle_type_type type) override;

    const Vec3 origin;
    const double radius, forceConstant;
};

class SphereOut : public PotentialOrder1 {
    using super = PotentialOrder1;
public:
    SphereOut(const std::string& particleType, double forceConstant, const Vec3& origin, double radius);

    virtual double getRelevantLengthScale() const noexcept override;

    virtual double getMaximalForce(double kbt) const noexcept override;

    double calculateEnergy(const Vec3 &position) const override;

    void calculateForce(Vec3 &force, const Vec3 &position) const override;

    void calculateForceAndEnergy(Vec3 &force, double &energy, const Vec3 &position) const override;

    void describe(std::ostream &os) const override;
protected:
    friend class readdy::model::KernelContext;

    void configureForType(const KernelContext *const ctx, const particle_type_type type) override;

    const Vec3 origin;
    const double radius, forceConstant;
};

template<typename T>
const std::string getPotentialName(typename std::enable_if<std::is_base_of<Cube, T>::value>::type * = 0) {
    return "Cube";
}
template<typename T>
const std::string getPotentialName(typename std::enable_if<std::is_base_of<SphereIn, T>::value>::type* = 0) {
    return "SphereIn";
}
template <typename T>
const std::string getPotentialName(typename std::enable_if<std::is_base_of<SphereOut, T>::value>::type* = 0) {
    return "SphereOut";
}

NAMESPACE_END(potentials)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
