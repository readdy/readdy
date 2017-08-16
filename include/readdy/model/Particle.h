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
#include "Vec3.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)

class READDY_API Particle {
public:

    using id_type = unsigned long;
    using pos_type = Vec3;
    using type_type = particle_type_type;

    Particle(scalar x, scalar y, scalar z, type_type type);

    Particle(Vec3 pos, type_type type);

    Particle(Vec3 pos, type_type type, id_type id);

    Particle(const Particle&) = default;

    Particle& operator=(const Particle&) = default;

    Particle(Particle&&) = default;

    Particle& operator=(Particle&&) = default;

    virtual ~Particle();

    const Vec3 &getPos() const;

    Vec3 &getPos();

    const type_type &getType() const;

    const id_type getId() const;

    bool operator==(const Particle &rhs) const;

    bool operator!=(const Particle &rhs) const;

    friend std::ostream &operator<<(std::ostream &, const Particle &);

    static id_type nextId();

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
