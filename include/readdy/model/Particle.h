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

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)

class READDY_API Particle {
public:

    using id_type = unsigned long;
    using pos_type = Vec3;
    using type_type = ParticleTypeId;

    Particle(scalar x, scalar y, scalar z, type_type type)
            : id(std::atomic_fetch_add<unsigned long>(&id_counter, 1L)), pos(x, y, z), type(type) {};

    Particle(Vec3 pos, type_type type) : pos(pos), type(type), id(std::atomic_fetch_add<id_type>(&id_counter, 1)) {}

    Particle(Vec3 pos, type_type type, id_type id) : pos(pos), type(type), id(id) {}

    Particle(const Particle&) = default;

    Particle& operator=(const Particle&) = default;

    Particle(Particle&&) = default;

    Particle& operator=(Particle&&) = default;

    virtual ~Particle() = default;

    const Vec3 &getPos() const {
        return pos;
    }

    Vec3 &getPos() {
        return pos;
    }

    const type_type &getType() const {
        return type;
    }

    const id_type getId() const {
        return id;
    }

    bool operator==(const Particle &rhs) const {
        return rhs.id == id;
    }

    bool operator!=(const Particle &rhs) const {
        return !(*this == rhs);
    }

    friend std::ostream &operator<<(std::ostream &os, const Particle &p) {
        os << "Particle(id=" << p.id << ", type=" << p.type << ", pos=" << p.pos << ")";
        return os;
    }

    static id_type nextId() {
        return std::atomic_fetch_add<id_type>(&id_counter, 1);
    }

protected:
    Vec3 pos;
    type_type type;
    id_type id;

    static std::atomic<id_type> id_counter;
};

class TopologyParticle : public Particle {
public:
    using Particle::Particle;
};

NAMESPACE_END(model)
NAMESPACE_END(readdy)
