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

#include <readdy/common/logging.h>
#include <readdy/kernel/singlecpu/model/SCPUParticleData.h>
#include <readdy/kernel/mpi/model/MPIDomain.h>

namespace readdy::kernel::mpi {

/**
 * An Entry similar to SCPU with additional MPI related data.
 * todo maybe ask particle instead of NL whether it is in halo (measure performance)
 */
struct MPIEntry {
    using Particle = readdy::model::Particle;
    using Force = Vec3;

    explicit MPIEntry(const Particle &particle, bool responsible = true, int rank = -1)
            : pos(particle.pos()), force(Force()), type(particle.type()),
              id(particle.id()), deactivated(false), responsible(responsible), rank(rank) {

    }

    [[nodiscard]] bool is_deactivated() const {
        return deactivated;
    }

    [[nodiscard]] const Particle::Position &position() const {
        return pos;
    }

    Particle::Position pos;
    Force force;
    ParticleId id;
    ParticleTypeId type;
    bool deactivated;
    /**
     * rank, which is the MPI rank this particle belongs to.
     * Usually is determined by position but not necessarily, it states which rank is responsible for this.
     * In evaluating forces and reactions it helps to know the rank without parsing the domain object (more deref).
     * Also the compound (rank, id) might provide a unique identifier if necessary (not planning to).
     * Helpful when avoiding double counting of pair-interaction observables e.g. virial across domain boundaries.
     */
    int rank;
    /**
     * The responsible flag indicates that this particle will be sent to other workers during sync,
     * and that it will not be deleted when applying received sync data from other workers.
     * States that the worker associated with this kernel/rank is responsible for this particle.
     * Particles, that the rank is not responsible for, are deleted/deactivated/replaced during sync.
     *
     * Responsibility depends on the operation to be performed:
     * - Initially responsibility is responsible=isInDomainCore(p.pos)
     * - For diffusion step: if the initial position of the particle was in the domain core,
     *   then the present worker is responsible
     * - For reactions within domain core (i.e. excluding bimolecular reactions with particles in halo):
     *   same as diffusion, any particle that was in domain core has responsible=true,
     *   additionally all newly created particles have responsible=true regardless of position
     * - For reactions across domains: all possible reaction events are bimolecular and between particles of
     *   different responsibility, the worker whose (responsible) particle p1 has lower rankOfPosition(p.pos)
     *   will become responsible for both particles, while the worker whose (responsible) particle p2 has higher rank
     *   will drop the responsibility for p2. Then the lower rank worker is responsible for sending the updated state
     *   and the higher rank will drop its own p2 during synchronization.
     */
    bool responsible;
};

class MPIParticleData : public readdy::kernel::scpu::model::SCPUParticleData<MPIEntry> {
public:
    explicit MPIParticleData(const readdy::kernel::mpi::model::MPIDomain *domain) : SCPUParticleData(), _domain(domain) {}

    // additionally sets the `rank` and `responsible` fields
    void addParticles(const std::vector<Particle> &particles) override {
        for(const auto& p : particles) {
            MPIEntry entry {p};
            entry.rank = _domain->rankOfPosition(entry.pos);
            entry.responsible = (entry.rank == _domain->rank());
            if(!_blanks.empty()) {
                const auto idx = _blanks.back();
                _blanks.pop_back();
                entries.at(idx) = entry;
            } else {
                entries.emplace_back(entry);
            }
        }
    }

private:
    const readdy::kernel::mpi::model::MPIDomain *_domain;
};

using MPIDataContainer = readdy::kernel::mpi::MPIParticleData;

}
