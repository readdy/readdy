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

namespace readdy::kernel::mpi::util {

struct ParticlePOD {
    Vec3 position;
    ParticleTypeId typeId;

    ParticlePOD() : position(Vec3()), typeId(0) {}
    ParticlePOD(Vec3 position, ParticleTypeId typeId) : position(position), typeId(typeId) {}
    explicit ParticlePOD(const MPIEntry& mpiEntry) : position(mpiEntry.pos), typeId(mpiEntry.type) {}
};

enum tags {
    transmitObjects
};

template <typename T>
inline std::vector<T> receiveObjects(int senderRank, const MPI_Comm &comm) {
    MPI_Status status;
    MPI_Probe(senderRank, tags::transmitObjects, comm, &status);
    int byteCount;
    MPI_Get_count(&status, MPI_BYTE, &byteCount);
    const int number = byteCount / sizeof(T);
    std::vector<T> objects(number);
    MPI_Recv((void *) objects.data(), byteCount, MPI_BYTE, senderRank, tags::transmitObjects, comm,
             MPI_STATUS_IGNORE);
    return objects;
}

template <typename T>
inline void receiveAppendObjects(int senderRank, std::vector<T> &result, const MPI_Comm &comm) {
    MPI_Status status;
    MPI_Probe(senderRank, tags::transmitObjects, comm, &status);
    int byteCount;
    MPI_Get_count(&status, MPI_BYTE, &byteCount);
    const int number = byteCount / sizeof(T);
    const std::size_t sizeBefore = result.size();
    result.resize(sizeBefore + number);
    MPI_Recv((void *) (result.data() + sizeBefore), byteCount, MPI_BYTE, senderRank, tags::transmitObjects, comm,
             MPI_STATUS_IGNORE);
}

template <typename T>
inline void sendObjects(int targetRank, std::vector<T> objects, const MPI_Comm &comm) {
    MPI_Send((void *) objects.data(), static_cast<int>(objects.size() * sizeof(T)), MPI_BYTE,
             targetRank, tags::transmitObjects, comm);
}

inline void sendThenReceive(std::array<std::size_t, 3> otherDirection, std::vector<util::ParticlePOD> &own,
                     std::vector<util::ParticlePOD> &other, const model::MPIDomain &domain, const MPI_Comm &comm) {
    const auto otherFlatIndex = domain.neighborIndex.index(otherDirection);
    const auto nType = domain.neighborTypes().at(otherFlatIndex);
    if (nType == model::MPIDomain::NeighborType::regular) {
        const auto otherRank = domain.neighborRanks().at(otherFlatIndex);
        // send
        std::vector<util::ParticlePOD> objects;
        objects.insert(objects.begin(), own.begin(), own.end());
        objects.insert(objects.begin(), other.begin(), other.end());
        util::sendObjects(otherRank, objects, comm);
        // receive, todo consider receiveAppendObjects
        auto received = util::receiveObjects<util::ParticlePOD>(otherRank, comm);
        other.insert(other.end(), received.begin(), received.end());
    }
}

inline void receiveThenSend(std::array<std::size_t, 3> otherDirection, std::vector<util::ParticlePOD> &own,
                     std::vector<util::ParticlePOD> &other, const model::MPIDomain &domain, const MPI_Comm &comm) {
    const auto otherFlatIndex = domain.neighborIndex.index(otherDirection);
    const auto nType = domain.neighborTypes().at(otherFlatIndex);
    if (nType == model::MPIDomain::NeighborType::regular) {
        const auto otherRank = domain.neighborRanks().at(otherFlatIndex);
        // receive
        auto received = util::receiveObjects<util::ParticlePOD>(otherRank, comm);
        // send
        std::vector<util::ParticlePOD> objects;
        objects.insert(objects.begin(), own.begin(), own.end());
        objects.insert(objects.begin(), other.begin(), other.end());
        util::sendObjects(otherRank, objects, comm);

        // append the received particles only after sending so we do not give the sender back its own particles
        other.insert(other.end(), received.begin(), received.end());
    }
}

}
