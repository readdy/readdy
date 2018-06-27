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
 * Declaration of the base class for all order 2 potentials. They basically have calculateForce, calculateEnergy and
 * calculateForceAndEnergy methods, which take a modifiable reference argument and the difference vector x_ij between
 * two particles.
 * Further, subclasses have to implement getCutoffRadius so that the neighbor list can be created more efficiently.
 *
 * @file PotentialOrder2.h
 * @brief Declaration of the base class for all order 2 potentials.
 * @author clonker
 * @date 31.05.16
 */

#pragma once

#include <cmath>

#include "Potential.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
class ParticleTypeRegistry;

NAMESPACE_BEGIN(potentials)
class PotentialRegistry;
class PotentialOrder2 : public Potential {
protected:
    using particle_type_type = readdy::model::Particle::type_type;
public:
    PotentialOrder2(particle_type_type type1, particle_type_type type2)
            : _particleType1(type1), _particleType2(type2) {}

    virtual scalar calculateEnergy(const Vec3 &x_ij) const = 0;

    virtual void calculateForce(Vec3 &force, const Vec3 &x_ij) const = 0;

    void calculateForceAndEnergy(Vec3 &force, scalar &energy, const Vec3 &x_ij) const {
        energy += calculateEnergy(x_ij);
        calculateForce(force, x_ij);
    };

    scalar getCutoffRadius() const {
        return std::sqrt(getCutoffRadiusSquared());
    };

    virtual scalar getCutoffRadiusSquared() const = 0;

    friend std::ostream &operator<<(std::ostream &os, const PotentialOrder2 &potential) {
        os << potential.describe();
        return os;
    }

    particle_type_type particleType1() const {
        return _particleType1;
    }

    particle_type_type particleType2() const {
        return _particleType2;
    }

protected:
    particle_type_type _particleType1, _particleType2;
};

NAMESPACE_END(potentials)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
