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
 * << detailed description >>
 *
 * @file TrajectoryEntry.h
 * @brief << brief description >>
 * @author clonker
 * @date 16.03.17
 * @copyright BSD-3
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
            : typeId(p.type()), id(p.id()), pos(p.pos()), flavor(ptr.infoOf(p.type()).flavor) {}

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
