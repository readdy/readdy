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
 * << detailed description >>
 *
 * @file TrajectoryEntry.h
 * @brief << brief description >>
 * @author clonker
 * @date 16.03.17
 * @copyright GNU Lesser General Public License v3.0
 */

#pragma once

#include <spdlog/fmt/ostr.h>
#include <readdy/common/common.h>
#include <readdy/model/Particle.h>
#include <readdy/model/ParticleTypeRegistry.h>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(observables)

struct TrajectoryEntry {

    TrajectoryEntry() = default;

    explicit TrajectoryEntry(const readdy::model::Particle &p, const readdy::model::ParticleTypeRegistry& ptr)
            : typeId(p.getType()), id(p.getId()), pos(p.getPos()), flavor(ptr.infoOf(p.getType()).flavor) {}

    readdy::model::Particle::type_type typeId {0};
    readdy::model::Particle::id_type id {0};
    readdy::model::particle_flavor flavor {0};
    readdy::model::Particle::pos_type pos;

    friend std::ostream &operator<<(std::ostream & /*os*/, const TrajectoryEntry & /*p*/);
};

inline std::ostream &operator<<(std::ostream &os, const TrajectoryEntry &p) {
    os << "TrajectoryEntry(id=" << p.id << ", type=" << p.typeId << ", position=" << p.pos << ", flavor=" << (int) p.flavor
       << ")";
    return os;
}

NAMESPACE_END(observables)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
