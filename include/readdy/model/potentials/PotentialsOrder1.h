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

#ifndef READDY_MAIN_POTENTIALSORDER1_H
#define READDY_MAIN_POTENTIALSORDER1_H

namespace readdy {
namespace model {

namespace potentials {

class CubePotential : public PotentialOrder1 {
    using super = PotentialOrder1;
public:
    CubePotential(const std::string& particleType, double forceConstant, const Vec3& origin, const Vec3& extent,
                  bool considerParticleRadius = true);

    const Vec3 &getOrigin() const;

    const Vec3 &getExtent() const;

    double getForceConstant() const;

    bool isConsiderParticleRadius() const;

    double getParticleRadius() const;

    virtual double getRelevantLengthScale() const noexcept override;

    virtual double getMaximalForce(double kbt) const noexcept override;

    std::string describe() override;

    friend std::ostream &operator<<(std::ostream &os, const CubePotential &potential);

    double calculateEnergy(const Vec3 &position) const override;

    void calculateForce(Vec3 &force, const Vec3 &position) const override;

    void calculateForceAndEnergy(Vec3 &force, double &energy, const Vec3 &position) const override;

protected:
    friend class readdy::model::KernelContext;

    void configureForType(const KernelContext *const ctx, const unsigned int type) override;

    const Vec3 origin, extent, min, max;
    const double forceConstant;
    const bool considerParticleRadius;
    double particleRadius;
};

class SpherePotential : public PotentialOrder1 {
    using super = PotentialOrder1;
public:
    SpherePotential(const std::string& particleType, double forceConstant, const Vec3& origin, double radius);

    const Vec3 &getOrigin() const;
    double getRadius() const;
    double getForceConstant() const;

    virtual double getRelevantLengthScale() const noexcept override;

    virtual double getMaximalForce(double kbt) const noexcept override;

    std::string describe() override;

    friend std::ostream &operator<<(std::ostream &os, const SpherePotential &potential);

    double calculateEnergy(const Vec3 &position) const override;

    void calculateForce(Vec3 &force, const Vec3 &position) const override;

    void calculateForceAndEnergy(Vec3 &force, double &energy, const Vec3 &position) const override;

protected:
    friend class readdy::model::KernelContext;

    void configureForType(const KernelContext *const ctx, const unsigned int type) override;

    const Vec3 origin;
    const double radius, forceConstant;
};

template<typename T>
const std::string getPotentialName(typename std::enable_if<std::is_base_of<CubePotential, T>::value>::type * = 0) {
    return "Cube";
}
template<typename T>
const std::string getPotentialName(typename std::enable_if<std::is_base_of<SpherePotential, T>::value>::type* = 0) {
    return "Sphere";
}
}
}
}

#endif //READDY_MAIN_POTENTIALSORDER1_H
