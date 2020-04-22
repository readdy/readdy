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
 * @file MPIStateModel.cpp
 * @brief « brief description »
 * @author chrisfroe
 * @date 28.05.19
 */

#include <readdy/kernel/mpi/MPIStateModel.h>
#include <readdy/common/Timer.h>

namespace readdy::kernel::mpi {

MPIStateModel::MPIStateModel(Data &data, const readdy::model::Context &context, const model::MPIDomain *domain)
        : _data(data), _context(context), _domain(domain), _neighborList(data, context, domain) {}

std::vector<readdy::Vec3> MPIStateModel::getParticlePositions() const {
    const auto data = getParticleData();
    std::vector<Vec3> target{};
    target.reserve(data->size());
    for (const auto &entry : *data) {
        if (!entry.deactivated) target.push_back(entry.pos);
    }
    return target;
}

// encapsulate the following combination of Gather and Gatherv, e.g. for gathering particles or observables
std::vector<MPIStateModel::Particle>
MPIStateModel::gatherParticles() const {
    if (_domain->isIdleRank()) {
        return {};
    }
    readdy::util::Timer timer("MPIStateModel::gatherParticles");
    auto &data = _data.get();

    // find out how many particles (and bytes) each worker sends
    int nParticles = 0;
    if (_domain->isWorkerRank()) {
        nParticles = std::count_if(data.begin(), data.end(), [](const MPIEntry &entry) {return not entry.deactivated and entry.responsible;});
    }
    std::vector<int> numberParticles(_domain->nUsedRanks(), 0);

    MPI_Gather(&nParticles, 1, MPI_INT, numberParticles.data(), 1, MPI_INT, 0, _commUsedRanks);

    std::vector<int> numberBytes(_domain->nUsedRanks(), 0);
    // convert number of particles to number of bytes from each rank
    for (std::size_t i = 0; i<numberParticles.size(); ++i) {
        numberBytes[i] = static_cast<int>(numberParticles[i] * sizeof(util::ParticlePOD));
    }

    // now send and receive thin particles
    if (_domain->isMasterRank()) {
        std::size_t totalNumberParticles = std::accumulate(numberParticles.begin(), numberParticles.end(), 0);
        std::vector<util::ParticlePOD> thinParticles(totalNumberParticles, {{0., 0., 0.}, 0}); // receive buffer
        std::vector<int> displacements(_domain->nUsedRanks(), 0);
        displacements[0] = 0;
        for (int i = 1; i<displacements.size(); ++i) {
            displacements[i] = displacements[i-1] + numberBytes[i-1];
        }
        MPI_Gatherv(nullptr, 0, MPI_BYTE, thinParticles.data(), numberBytes.data(), displacements.data(), MPI_BYTE, 0, _commUsedRanks);
        // convert to particles
        std::vector<Particle> particles;
        std::for_each(thinParticles.begin(), thinParticles.end(),
                      [&particles](const util::ParticlePOD &tp) {
                          particles.emplace_back(tp.position, tp.typeId);
                      });
        return particles;
    } else {
        // prepare send data
        std::vector<util::ParticlePOD> thinParticles;
        for (const MPIEntry &entry : _data.get()) {
            if (not entry.deactivated and entry.responsible) {
                thinParticles.emplace_back(entry);
            }
        }
        MPI_Gatherv((void *) thinParticles.data(), static_cast<int>(thinParticles.size() * sizeof(util::ParticlePOD)), MPI_BYTE, nullptr, nullptr, nullptr, nullptr, 0, _commUsedRanks);
        return {};
    }
}

void MPIStateModel::resetReactionCounts() {
    if (!reactionCounts().empty()) {
        for (auto &e : reactionCounts()) {
            e.second = 0;
        }
    } else {
        const auto &reactions = _context.get().reactions();
        for (const auto &entry : reactions.order1()) {
            for (auto reaction : entry.second) {
                reactionCounts()[reaction->id()] = 0;
            }
        }
        for (const auto &entry : reactions.order2()) {
            for (auto reaction : entry.second) {
                reactionCounts()[reaction->id()] = 0;
            }
        }
    }
}

void MPIStateModel::toDenseParticleIndices(
        std::vector<std::size_t>::iterator begin,
        std::vector<std::size_t>::iterator end) const {
    const auto &blanks = _data.get().blanks();
    std::transform(begin, end, begin, [&blanks](const std::size_t &ix) {
        auto result = ix;
        for (auto blankIx : blanks) {
            if (blankIx < ix) {
                --result;
            }
        }
        return result;
    });
}

void MPIStateModel::clear() {
    getParticleData()->clear();
    //topologies().clear();
    reactionRecords().clear();
    resetReactionCounts();
    virial() = {};
    energy() = 0;
}

void MPIStateModel::distributeParticle(const Particle &p) {
    distributeParticles({p});
}

void MPIStateModel::addParticles(const std::vector<Particle> &particles) {
    getParticleData()->addParticles(particles);
}

// todo use MPI_Scatter
void MPIStateModel::distributeParticles(const std::vector<Particle> &ps) {
    if (_domain->isIdleRank()) {
        return;
    }
    readdy::util::Timer timer("MPIStateModel::distributeParticles");
    if (_domain->isMasterRank()) {
        std::unordered_map<int, std::vector<util::ParticlePOD>> targetParticleMap;
        for (const auto &particle : ps) {
            int target = _domain->rankOfPosition(particle.pos());
            assert(target < domain()->nUsedRanks());
            assert(target != 0);
            const auto &find = targetParticleMap.find(target);
            if (find != targetParticleMap.end()) {
                find->second.emplace_back(particle.pos(), particle.type());
            } else {
                targetParticleMap.emplace(std::make_pair(
                        target, std::vector<util::ParticlePOD>{{particle.pos(), particle.type()}}
                ));
            }
        }
        std::vector<char> isTarget(_domain->nUsedRanks(), false);
        for (auto&&[target, vec] : targetParticleMap) {
            isTarget[target] = true;
        }
        // broadcast target
        MPI_Bcast(isTarget.data(), isTarget.size(), MPI_CHAR, 0, _commUsedRanks);
        // send
        for (auto&&[target, thinParticles] : targetParticleMap) {
            MPI_Request req;
            MPI_Isend((void *) thinParticles.data(),
                      static_cast<int>(thinParticles.size() * sizeof(util::ParticlePOD)), MPI_BYTE,
                      target, util::tags::transmitObjects, _commUsedRanks, &req);
            MPI_Request_free(&req);
        }
        MPI_Barrier(_commUsedRanks);
    } else {
        std::vector<char> isTarget(_domain->nUsedRanks(), false);
        // am i one of the targets?
        MPI_Bcast(isTarget.data(), isTarget.size(), MPI_CHAR, 0, _commUsedRanks);
        if (isTarget[_domain->rank()]) {
            // receive, ignore argument p here
            const auto thinParticles = util::receiveObjects<util::ParticlePOD>(0, _commUsedRanks);
            std::vector<Particle> particles;
            std::for_each(thinParticles.begin(), thinParticles.end(),
                          [&particles](const util::ParticlePOD &tp) {
                              particles.emplace_back(tp.position, tp.typeId);
                          });
            addParticles(particles);
        }
        MPI_Barrier(_commUsedRanks);
    }
}

std::vector<readdy::model::Particle> MPIStateModel::getParticles() const {
    const auto *data = getParticleData();
    std::vector<readdy::model::Particle> result;
    result.reserve(data->size());
    for (const auto &entry : *data) {
        if (!entry.deactivated) {
            result.push_back(data->toParticle(entry));
        }
    }
    return result;
}

void MPIStateModel::addParticle(const MPIStateModel::Particle &p) {
    addParticles({p});
}

// todo use mpi built in cartesian graph communicator and neighborhood collectives
// MPI_Neighbor_allgather(const void* sendbuf, int sendcount,
//                        MPI_Datatype sendtype, void* recvbuf, int recvcount,
//                        MPI_Datatype recvtype, MPI_Comm comm)
void MPIStateModel::synchronizeWithNeighbors() {
    if (domain()->isIdleRank() or domain()->isMasterRank()) {
        return;
    }
    readdy::util::Timer timer("MPIStateModel::synchronizeWithNeighbors");
    auto& data = _data.get();
    std::vector<util::ParticlePOD> own; // particles that this worker is responsible for
    std::vector<std::size_t> removedEntries; // particles that this worker is NOT responsible for

    // gather own responsible and prepare data structure
    // i.e. gather to-be-removed indices,
    // and re-tag particles that are currently responsible but not in core of domain
    for (size_t i = 0; i < data.size(); ++i) {
        MPIEntry& entry = data.entry_at(i);
        if (not entry.deactivated and entry.responsible) {
            own.emplace_back(entry);
            if (domain()->isInDomainHalo(entry.pos)) {
                entry.responsible = false;
            }
        } else if (not entry.deactivated and not entry.responsible) {
            removedEntries.push_back(i);
        }
    }

    readdy::util::Timer t1("MPIStateModel::synchronizeWithNeighbors.plimpton");
    const auto &pbc = _context.get().periodicBoundaryConditions();
    // Plimpton synchronization
    std::vector<util::ParticlePOD> other; // particles received by other workers
    for (unsigned int coord=0; coord<3; coord++) { // east-west, north-south, up-down
        const auto idx = domain()->myIdx()[coord];
        if (idx % 2 == 0) {
            // send + then receive +
            std::array<std::size_t, 3> otherDirection {1,1,1}; // (1,1,1) is self
            otherDirection.at(coord) += 1;
            auto received1 = util::sendThenReceive(otherDirection, own, other, *domain(), commUsedRanks());

            // receive - then send -
            otherDirection = {1,1,1};
            otherDirection.at(coord) -= 1;
            std::vector<util::ParticlePOD> received2;
            if (domain()->nDomainsPerAxis()[coord] == 2 and pbc[coord]) {
                // skip because we have already communicated with that one
            } else {
                received2 = util::receiveThenSend(otherDirection, own, other, *domain(), commUsedRanks());
            }
            // after data from both directions have been received we can merge them with `other`,
            // so they will be communicated along other coordinates
            other.insert(other.end(), received1.begin(), received1.end());
            other.insert(other.end(), received2.begin(), received2.end());
        } else {
            // receive - then send -
            std::array<std::size_t, 3> otherDirection {1,1,1};
            otherDirection.at(coord) -= 1;
            auto received1 = util::receiveThenSend(otherDirection, own, other, *domain(), commUsedRanks());

            // send + then receive +
            otherDirection = {1,1,1};
            otherDirection.at(coord) += 1;
            std::vector<util::ParticlePOD> received2;
            if (domain()->nDomainsPerAxis()[coord] == 2 and pbc[coord]) {
                // skip because we have already communicated with that one
            } else {
                received2 = util::sendThenReceive(otherDirection, own, other, *domain(), commUsedRanks());
            }

            other.insert(other.end(), received1.begin(), received1.end());
            other.insert(other.end(), received2.begin(), received2.end());
        }
    }
    t1.stop();

    // only add new entries if in domain coreOrHalo and additionally set responsible=true if in core
    std::vector<MPIEntry> newEntries;
    for (const auto &p : other) {
        if (domain()->isInDomainCore(p.position)) {
            // gets added and worker is responsible
            Particle particle(p.position, p.typeId);
            MPIEntry entry(particle, true, domain()->rank());
            newEntries.emplace_back(entry);
        } else if (domain()->isInDomainCoreOrHalo(p.position)) {
            // gets added but worker is not responsible
            Particle particle(p.position, p.typeId);
            MPIEntry entry(particle, false, domain()->rankOfPosition(p.position));
            newEntries.emplace_back(entry);
        } else {
            // does not get added
        }
    }
    auto update = std::make_pair(std::move(newEntries), std::move(removedEntries));
    data.update(std::move(update));
}

}
