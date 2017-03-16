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
#include <ostream>
#include "Vec3.h"

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)

class READDY_API Particle {
public:

    using id_type = unsigned long;
    using pos_type = Vec3;
    using type_type = unsigned short;
    using flavor_t = std::uint8_t;

    static constexpr flavor_t FLAVOR_NORMAL = 0;
    static constexpr flavor_t FLAVOR_TOPOLOGY = 1;
    static constexpr flavor_t FLAVOR_MEMBRANE = 2;

    Particle(double x, double y, double z, type_type type);

    Particle(Vec3 pos, type_type type);

    Particle(Vec3 pos, type_type type, id_type id);

    virtual ~Particle();

    const Vec3 &getPos() const;

    Vec3 &getPos();

    flavor_t getFlavor() const;

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
    flavor_t flavor;

    static std::atomic<id_type> id_counter;
};

class TopologyParticle : public Particle {
    using super = Particle;
public:
    TopologyParticle(double x, double y, double z, type_type type);

    TopologyParticle(Vec3 pos, type_type type);

    TopologyParticle(Vec3 pos, type_type type, id_type id);

};

NAMESPACE_END(model)
NAMESPACE_END(readdy)
