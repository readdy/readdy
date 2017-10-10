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
 * Declaration of the base class of all order 1 potentials.
 * Subclasses have to implement calculateEnergy, calculateForce and calculateForceAndEnergy.
 * The first three methods take a modifiable reference and a particle's position. The last method is for replication
 * of the potential, so that it can be assigned to multiple particle types.
 *
 * @file PotentialOrder1.h
 * @brief Declaration of the base class of all order 1 potentials
 * @author clonker
 * @date 31.05.16
 */

#pragma once
#include <readdy/model/potentials/Potential.h>

#include <utility>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
class ParticleTypeRegistry;

NAMESPACE_BEGIN(potentials)
class PotentialRegistry;
class PotentialOrder1 : public Potential {
protected:
    using particle_type_type = readdy::model::Particle::type_type;
public:
    explicit PotentialOrder1(particle_type_type ptype) : Potential(1), _particleType(ptype) {}

    virtual scalar calculateEnergy(const Vec3 &position) const = 0;

    virtual void calculateForce(Vec3 &force, const Vec3 &position) const = 0;

    virtual void calculateForceAndEnergy(Vec3 &force, scalar &energy, const Vec3 &position) const = 0;

    virtual scalar getRelevantLengthScale() const noexcept = 0;

    friend std::ostream &operator<<(std::ostream &os, const PotentialOrder1 &potential) {
        os << potential.describe();
        return os;
    }

    particle_type_type particleType() const {
        return _particleType;
    }

protected:
    particle_type_type _particleType;
};

NAMESPACE_END(potentials)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
