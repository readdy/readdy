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
    explicit PotentialOrder1(particle_type_type ptype) : _particleType(ptype) {}

    virtual scalar calculateEnergy(const Vec3 &position) const = 0;

    virtual void calculateForce(Vec3 &force, const Vec3 &position) const = 0;

    void calculateForceAndEnergy(Vec3 &force, scalar &energy, const Vec3 &position) const {
        energy += calculateEnergy(position);
        calculateForce(force, position);
    };

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
