/********************************************************************
 * Copyright © 2019 Computational Molecular Biology Group,          *
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
 * @file MPIParticleData.h
 * @brief Header file defining the type of entry, which the data container uses.
 * @author chrisfroe
 * @date 13.06.19
 */

#pragma once

#include <readdy/kernel/singlecpu/model/SCPUParticleData.h>

namespace readdy::kernel::mpi {

/**
 * An Entry similar to SCPU with additional rank, which is the MPI rank this particle belongs to, 
 * because IDs are local to kernel and we need a unique identifier for particles, 
 * which is the compound (rank, id).
 */
struct MPIEntry {
    using Particle = readdy::model::Particle;
    using Force = Particle::Position;

    explicit MPIEntry(const Particle &particle, int rank = -1)
            : pos(particle.pos()), force(Force()), type(particle.type()), deactivated(false),
              id(particle.id()), rank(rank) {}

    bool is_deactivated() const {
        return deactivated;
    }

    const Particle::Position &position() const {
        return pos;
    }

    Force force;
    Particle::Position pos;
    Particle::Id id;
    Particle::TypeId type;
    bool deactivated;
    int rank;
};

using MPIDataContainer = readdy::kernel::scpu::model::SCPUParticleData<MPIEntry>;

}
