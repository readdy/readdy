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

namespace readdy::kernel::mpi {

MPIStateModel::MPIStateModel(Data &data, const readdy::model::Context &context) : _data(data), _context(context) {}

const std::vector<readdy::Vec3> MPIStateModel::getParticlePositions() const {
    throw std::runtime_error("impl");
}

// todo procedure to send variable type objects?
const std::vector<MPIStateModel::Particle>
MPIStateModel::getParticles() const {
    if (not mpiutil::isRequiredRank(*_domain)) {
        return {};
    }
    // todo test this
    // todo time this
    // find out how many particles (and bytes) each worker sends
    int nParticles = 0;
    if (_domain->rank != 0) {
        nParticles = _data.get().size();
    }
    std::vector<int> numberBytes(_domain->worldSize, 0);
    MPI_Gather(&nParticles, 1, MPI_INT, &numberBytes, numberBytes.size(), MPI_INT, 0, MPI_COMM_WORLD);
    // convert number of particles to number of bytes from each rank
    for (auto &n : numberBytes) {
        n = static_cast<int>(n * sizeof(mpiutil::ThinParticle));
    }

    // now send and receive thin particles
    if (_domain->rank == 0) {
        std::vector<mpiutil::ThinParticle> thinParticles; // receive buffer
        std::vector<int> displacements(_domain->worldSize, 0);
        displacements[0] = 0;
        for (int i = 1; i<displacements.size(); ++i) {
            displacements[i] = displacements[i-1] + numberBytes[i-1];
        }
        MPI_Gatherv(nullptr, 0, MPI_BYTE, thinParticles.data(), numberBytes.data(), displacements.data(), MPI_BYTE, 0, MPI_COMM_WORLD);
        // convert to particles
        std::vector<Particle> particles;
        std::for_each(thinParticles.begin(), thinParticles.end(),
                      [&particles](const mpiutil::ThinParticle &tp) {
                          particles.emplace_back(tp.position, tp.typeId);
                      });
        return particles;
    } else {
        // prepare send data
        std::vector<mpiutil::ThinParticle> thinParticles;
        for (const auto &entry : _data.get()) {
            thinParticles.emplace_back(entry.position(), entry.type);
        }
        MPI_Gatherv((void *) thinParticles.data(), static_cast<int>(thinParticles.size() * sizeof(mpiutil::ThinParticle)), MPI_BYTE, nullptr, nullptr, nullptr, nullptr, 0, MPI_COMM_WORLD);
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

void MPIStateModel::addParticle(const Particle &p) {
    if (not mpiutil::isRequiredRank(*_domain)) {
        return;
    }
    readdy::log::trace("MPIStateModel::addParticle");
    MPI_Barrier(*_commUsedRanks);
    if (_domain->rank == 0) {
        int targetRank = _domain->rankOfPosition(p.pos());
        // broadcast target
        MPI_Bcast(&targetRank, 1, MPI_INT, 0, *_commUsedRanks);
        // send
        std::vector<mpiutil::ThinParticle> thinParticles{{p.pos(), p.type()}};
        MPI_Send((void *) thinParticles.data(),
                 static_cast<int>(thinParticles.size() * sizeof(mpiutil::ThinParticle)), MPI_BYTE,
                 targetRank, mpiutil::tags::sendParticles, *_commUsedRanks);
    } else {
        // am i the target?
        int targetRank;
        MPI_Bcast(&targetRank, 1, MPI_INT, 0, *_commUsedRanks);
        if (_domain->rank == targetRank) {
            // receive, ignore argument p here
            const auto thinParticles = mpiutil::receiveParticlesFrom(0, *_commUsedRanks);
            std::vector<Particle> particles;
            std::for_each(thinParticles.begin(), thinParticles.end(),
                          [&particles](const mpiutil::ThinParticle &tp) {
                              particles.emplace_back(tp.position, tp.typeId);
                          });
            getParticleData()->addParticles(particles);
        }
    }
};

void MPIStateModel::addParticles(const std::vector<Particle> &p) {
    if (not mpiutil::isRequiredRank(*_domain)) {
        return;
    }
    // todo use MPI_Scatter
    if (_domain->rank == 0) {
        std::unordered_map<int, std::vector<mpiutil::ThinParticle>> targetParticleMap;
        for (const auto &particle : p) {
            int target = _domain->rankOfPosition(particle.pos());
            const auto &find = targetParticleMap.find(target);
            if (find != targetParticleMap.end()) {
                find->second.emplace_back(particle.pos(), particle.type());
            } else {
                targetParticleMap.emplace(std::make_pair(
                        target, std::vector<mpiutil::ThinParticle>{{particle.pos(), particle.type()}}
                ));
            }

        }
        std::vector<char> isTarget(_domain->worldSize, false);
        for (auto&&[target, vec] : targetParticleMap) {
            isTarget[target] = true;
        }
        // broadcast target
        MPI_Bcast(&isTarget, isTarget.size(), MPI_CHAR, 0, *_commUsedRanks);
        // send
        for (auto&&[target, thinParticles] : targetParticleMap) {
            MPI_Send((void *) thinParticles.data(),
                     static_cast<int>(thinParticles.size() * sizeof(mpiutil::ThinParticle)), MPI_BYTE,
                     target, mpiutil::tags::sendParticles, *_commUsedRanks);
        }
    } else {
        std::vector<char> isTarget(_domain->worldSize, false);
        // am i one of the targets?
        MPI_Bcast(&isTarget, isTarget.size(), MPI_CHAR, 0, *_commUsedRanks);
        if (isTarget[_domain->rank]) {
            // receive, ignore argument p here
            const auto thinParticles = mpiutil::receiveParticlesFrom(0, *_commUsedRanks);
            std::vector<Particle> particles;
            std::for_each(thinParticles.begin(), thinParticles.end(),
                          [&particles](const mpiutil::ThinParticle &tp) {
                              particles.emplace_back(tp.position, tp.typeId);
                          });
            getParticleData()->addParticles(particles);
        }
    }
};

}
