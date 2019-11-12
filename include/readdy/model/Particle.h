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
 * A particle is a composite of position, type, and id. The type is an "unsigned int" value, which is mapped by
 * the KernelContext.
 *
 * @file Particle.h
 * @brief Header file containing the definitions for the particle class.
 * @author clonker
 * @date 19.04.16
 */

#pragma once
#include <array>
#include <string>
#include <atomic>
#include <spdlog/fmt/ostr.h>

#include <readdy/common/common.h>
#include "readdy/common/ReaDDyVec3.h"

namespace readdy::model {

class Particle {
public:

    using Position = Vec3;

    Particle(scalar x, scalar y, scalar z, ParticleTypeId type)
            : _id(std::atomic_fetch_add<unsigned long>(&idCounter, 1L)), _pos(x, y, z), _type(type) {};

    Particle(Vec3 pos, ParticleTypeId type) : _pos(pos), _type(type), _id(std::atomic_fetch_add<ParticleId>(&idCounter, 1)) {}

    Particle(Vec3 pos, ParticleTypeId type, ParticleId id) : _pos(pos), _type(type), _id(id) {}

    Particle(const Particle &) = default;

    Particle &operator=(const Particle &) = default;

    Particle(Particle &&) = default;

    Particle &operator=(Particle &&) = default;

    virtual ~Particle() = default;

    [[nodiscard]] const Vec3 &pos() const {
        return _pos;
    }

    Vec3 &pos() {
        return _pos;
    }

    [[nodiscard]] ParticleTypeId type() const {
        return _type;
    }

    [[nodiscard]] ParticleId id() const {
        return _id;
    }

    bool operator==(const Particle &rhs) const {
        return rhs._id == _id;
    }

    bool operator!=(const Particle &rhs) const {
        return !(*this == rhs);
    }

    friend std::ostream &operator<<(std::ostream &os, const Particle &p) {
        os << "Particle(id=" << p._id << ", type=" << p._type << ", pos=" << p._pos << ")";
        return os;
    }

    static ParticleId nextId() {
        return std::atomic_fetch_add<ParticleId>(&idCounter, 1);
    }

protected:
    Vec3 _pos;
    ParticleTypeId _type;
    ParticleId _id;

    static std::atomic<ParticleId> idCounter;
};

class TopologyParticle : public Particle {
public:
    using Particle::Particle;
};

}
