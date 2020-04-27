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
#include <readdy/common/Timer.h>

namespace readdy::kernel::mpi::util {

struct ParticlePOD {
    Vec3 position;
    ParticleTypeId typeId;

    ParticlePOD() : position(Vec3()), typeId(0) {}

    ParticlePOD(Vec3 position, ParticleTypeId typeId) : position(position), typeId(typeId) {}

    explicit ParticlePOD(const MPIEntry &mpiEntry) : position(mpiEntry.pos), typeId(mpiEntry.type) {}
    explicit ParticlePOD(const readdy::model::Particle &particle) : position(particle.pos()), typeId(particle.type()) {}

    bool operator==(const ParticlePOD& other) const {
        return (this->position == other.position) and (this->typeId == other.typeId);
    }
};

struct HashPOD {
    std::size_t operator()(const readdy::kernel::mpi::util::ParticlePOD &pod) const {
        std::size_t seed{0};
        readdy::util::hash::combine(seed, std::hash<readdy::scalar>{}(pod.position.x));
        readdy::util::hash::combine(seed, std::hash<readdy::scalar>{}(pod.position.y));
        readdy::util::hash::combine(seed, std::hash<readdy::scalar>{}(pod.position.z));
        readdy::util::hash::combine(seed, std::hash<ParticleTypeId>{}(pod.typeId));
        return seed;
    }
};

enum tags {
    transmitObjects
};

template<typename T>
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

// todo here flatbuffers could be useful
template<typename T>
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

template<typename T>
inline void sendObjects(int targetRank, const std::vector<T> &objects, const MPI_Comm &comm) {
    MPI_Send((void *) objects.data(), static_cast<int>(objects.size() * sizeof(T)), MPI_BYTE,
             targetRank, tags::transmitObjects, comm);
}

inline std::ostream &operator<<(std::ostream& os, readdy::kernel::mpi::model::MPIDomain::NeighborType n) {
    switch(n) {
        case readdy::kernel::mpi::model::MPIDomain::NeighborType::self: os << "self"; break;
        case readdy::kernel::mpi::model::MPIDomain::NeighborType::nan: os << "nan"; break;
        case readdy::kernel::mpi::model::MPIDomain::NeighborType::regular: os << "regular"; break;
    }
    return os;
}

inline std::vector<util::ParticlePOD>
sendThenReceive(std::array<std::size_t, 3> otherDirection, std::vector<util::ParticlePOD> &own,
                std::vector<util::ParticlePOD> &other, const model::MPIDomain &domain, const MPI_Comm &comm) {
    const auto otherFlatIndex = domain.neighborIndex.index(otherDirection);
    const auto nType = domain.neighborTypes().at(otherFlatIndex);
    const auto otherRank = domain.neighborRanks().at(otherFlatIndex);
    //readdy::log::trace("rank={}, sendThenReceive, otherRank {}, otherType {}", domain.rank(), otherRank, nType);
    if (nType == model::MPIDomain::NeighborType::regular) {
        // send
        std::vector<util::ParticlePOD> objects;
        objects.insert(objects.end(), own.begin(), own.end());
        objects.insert(objects.end(), other.begin(), other.end());
        readdy::util::Timer t1("sendThenReceive.sendObjects");
        util::sendObjects(otherRank, objects, comm);
        t1.stop();
        //readdy::log::trace("rank={}, sendThenReceive.sent", domain.rank());
        // receive
        readdy::util::Timer t2("sendThenReceive.receiveObjects");
        auto received = util::receiveObjects<util::ParticlePOD>(otherRank, comm);
        t2.stop();
        //readdy::log::trace("rank={}, sendThenReceive.received", domain.rank());
        return received;
    } else {
        return {};
    }
}

inline std::vector<util::ParticlePOD>
receiveThenSend(std::array<std::size_t, 3> otherDirection, std::vector<util::ParticlePOD> &own,
                std::vector<util::ParticlePOD> &other, const model::MPIDomain &domain, const MPI_Comm &comm) {
    const auto otherFlatIndex = domain.neighborIndex.index(otherDirection);
    const auto nType = domain.neighborTypes().at(otherFlatIndex);
    const auto otherRank = domain.neighborRanks().at(otherFlatIndex);
    //readdy::log::trace("rank={}, receiveThenSend, otherRank {}, otherType {}", domain.rank(), otherRank, nType);
    if (nType == model::MPIDomain::NeighborType::regular) {
        // receive
        readdy::util::Timer t1("receiveThenSend.receiveObjects");
        auto received = util::receiveObjects<util::ParticlePOD>(otherRank, comm);
        t1.stop();
        //readdy::log::trace("rank={}, receiveThenSend.received", domain.rank());
        // send
        std::vector<util::ParticlePOD> objects;
        objects.insert(objects.end(), own.begin(), own.end());
        objects.insert(objects.end(), other.begin(), other.end());
        readdy::util::Timer t2("receiveThenSend.sendObjects");
        util::sendObjects(otherRank, objects, comm);
        t2.stop();
        //readdy::log::trace("rank={}, receiveThenSend.sent", domain.rank());
        return received;
    } else {
        return {};
    }
}

// specialized version for MPI, todo remove topology?
template<typename ParticleContainer, typename EvaluateOnParticle, typename InteractionContainer,
        typename EvaluateOnInteraction, typename TopologyContainer, typename EvaluateOnTopology>
inline void evaluateOnContainers(ParticleContainer &&particleContainer,
                                 const EvaluateOnParticle &evaluateOnParticle,
                                 InteractionContainer &&interactionContainer,
                                 const EvaluateOnInteraction &evaluateOnInteraction,
                                 TopologyContainer &&topologyContainer,
                                 const EvaluateOnTopology &evaluateOnTopology) {
    // Evaluate on particles
    {
        std::for_each(particleContainer.begin(), particleContainer.end(), [&](auto &&entry){
            if (!entry.deactivated) {
                evaluateOnParticle(entry);
            }
        });
    }

    // Evaluate on interactions
    {
        interactionContainer.forAllPairs(evaluateOnInteraction);
    }

    // Evaluate on topologies
    {
        for (auto &&topology : topologyContainer) {
            if (!topology->isDeactivated()) {
                evaluateOnTopology(topology);
            }
        }
    }
}

/**
 * Wrapper around two calls Gather and Gatherv,
 * to find out how many objects each one sends (1),
 * then set the appropriate displacements (2),
 * and then gather all variable length data (3).
 *
 * @tparam T, the type of sent objects
 * @param objects, the vector of sent objects
 * @param root, the rank of the worker which will end up with the union of all sent objects
 * @param domain, domain object with rank information of current worker
 * @param comm, communicator for the set of workers
 * @return the union of all sent objects if on root worker, else an empty vector
 */
template<typename T>
inline std::vector<T> gatherObjects(const std::vector<T> &objects, int root, const model::MPIDomain &domain, const MPI_Comm &comm) {
    /// (1) find out how many each one sends. In principle the root can also send objects.
    int number = objects.size();
    std::vector<int> nPerRank(domain.nUsedRanks(), 0);
    MPI_Gather(&number, 1, MPI_INT, nPerRank.data(), 1, MPI_INT, root, comm);

    std::vector<int> nPerRankBytes;
    std::vector<int> displacements;
    std::vector<T> results;

    if (domain.rank() == root) {
        /// (2) find out number of bytes and displacements
        nPerRankBytes.resize(domain.nUsedRanks());
        // convert number of objects to number of bytes from each rank
        for (std::size_t i = 0; i<nPerRank.size(); ++i) {
            nPerRankBytes[i] = static_cast<int>(nPerRank[i] * sizeof(T));
        }

        // prepare receive buffer and displacements
        std::size_t totalNumber = std::accumulate(nPerRank.begin(), nPerRank.end(), 0);
        results.resize(totalNumber);
        displacements.resize(domain.nUsedRanks());
        displacements[0] = 0;
        for (int i = 1; i<displacements.size(); ++i) {
            displacements[i] = displacements[i-1] + nPerRankBytes[i-1];
        }
    }

    /// (3) gather results
    MPI_Gatherv((void *) objects.data(), static_cast<int>(objects.size() * sizeof(T)), MPI_BYTE, results.data(),
                nPerRankBytes.data(), displacements.data(), MPI_BYTE, root, comm);
    return results;
}

}
