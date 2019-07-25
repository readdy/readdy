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
 * « detailed description »
 *
 * @file MPIUtils.h
 * @brief « brief description »
 * @author chrisfroe
 * @date 24.07.19
 */

#pragma once

#include <string>
#include <mpi.h>
#include <vector>

namespace readdy::kernel::mpi::mpiutil {

struct ThinParticle {
    ThinParticle(const Vec3 &position, ParticleTypeId typeId) : position(position), typeId(typeId) {};
    Vec3 position;
    ParticleTypeId typeId;
};

enum tags {
    sendParticles
};

inline std::vector<ThinParticle> receiveParticlesFrom(int sender, const MPI_Comm &comm) {
    MPI_Status status;
    MPI_Probe(sender, tags::sendParticles, comm, &status);
    int byteCount;
    MPI_Get_count(&status, MPI_BYTE, &byteCount);
    const int nParticles = byteCount / sizeof(ThinParticle);
    std::vector<ThinParticle> particles(nParticles, {{0.,0.,0.}, 0});
    MPI_Recv((void *) particles.data(), byteCount, MPI_BYTE, sender, tags::sendParticles, comm,
             MPI_STATUS_IGNORE);
    if (nParticles > 0) {
        return particles;
    } else {
        throw std::runtime_error("y u no send particles?");
    }
}

inline bool isRequiredRank(const model::MPIDomain &domain) {
    return domain.amINeeded();
}

}
