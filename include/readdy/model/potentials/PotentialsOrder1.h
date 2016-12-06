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

#include "PotentialOrder1.h"

#ifndef READDY_MAIN_POTENTIALSORDER1_H
#define READDY_MAIN_POTENTIALSORDER1_H

namespace readdy {
namespace model {

class Kernel;
namespace potentials {

class CubePotential : public PotentialOrder1 {

public:
    CubePotential(const Kernel *const kernel);

    virtual void configureForType(const unsigned int type) override;

    const Vec3 &getOrigin() const;

    void setOrigin(const Vec3 &origin);

    const Vec3 &getExtent() const;

    void setExtent(const Vec3 &extent);

    double getForceConstant() const;

    void setForceConstant(double forceConstant);

    bool isConsiderParticleRadius() const;

    void setConsiderParticleRadius(bool considerParticleRadius);

    double getParticleRadius() const;

    virtual double getRelevantLengthScale() const noexcept override;

    virtual double getMaximalForce(double kbt) const noexcept override;

protected:
    const Kernel *const kernel;
    Vec3 origin;
    Vec3 extent;
    Vec3 min{0, 0, 0}, max{0, 0, 0};
    double particleRadius = 0;
    double forceConstant = 1;
    bool considerParticleRadius = true;
};

class SpherePotential : public PotentialOrder1 {
public:
    SpherePotential(const Kernel * kernel);

    virtual void configureForType(const unsigned int type) override;

    const Vec3 &getOrigin() const;
    void setOrigin(const Vec3 &origin);
    double getRadius() const;
    void setRadius(double radius);
    double getForceConstant() const;
    void setForceConstant(double forceConstant);

    virtual double getRelevantLengthScale() const noexcept override;

    virtual double getMaximalForce(double kbt) const noexcept override;

protected:
    const Kernel * const kernel;
    Vec3 origin;
    double radius = 1;
    double forceConstant = 1;
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
