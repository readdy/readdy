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
 * @file MPICalculateForces.cpp
 * @brief « brief description »
 * @author chrisfroe
 * @date 03.06.19
 */

#include <readdy/kernel/mpi/actions/MPIActions.h>
#include <readdy/common/algorithm.h>
#include <readdy/common/boundary_condition_operations.h>

namespace readdy::kernel::mpi::actions {

namespace detail {
template<bool COMPUTE_VIRIAL>
void computeVirial(const Vec3& r_ij, const Vec3 &force, Matrix33 &virial);

template<>
void computeVirial<true>(const Vec3& r_ij, const Vec3 &force, Matrix33 &virial) {
    virial += math::outerProduct<Matrix33>(-1.*r_ij, force);
}

template<>
void computeVirial<false>(const Vec3& /*r_ij*/, const Vec3 &/*force*/, Matrix33 &/*virial*/) {}
}

template<bool COMPUTE_VIRIAL>
void MPICalculateForces::performImpl() {
    const auto &context = kernel->context();
    auto &stateModel = kernel->getMPIKernelStateModel();
    auto &data = *stateModel.getParticleData();
    auto &neighborList = *stateModel.getNeighborList();

    stateModel.energy() = 0;
    stateModel.virial() = Matrix33{{{0, 0, 0, 0, 0, 0, 0, 0, 0}}};

    const auto &potentials = context.potentials();

    if (!potentials.potentialsOrder1().empty() || !potentials.potentialsOrder2().empty()) {
        std::for_each(data.begin(), data.end(), [](auto &entry) {
            entry.force = {0, 0, 0};
        });
    }

    auto order1eval = [&](auto &entry){
        for (const auto &po1 : potentials.potentialsOf(entry.type)) {
            po1->calculateForceAndEnergy(entry.force, stateModel.energy(), entry.position());
        }
    };

    // order 2 eval
    const auto &box = context.boxSize().data();
    const auto &pbc = context.periodicBoundaryConditions().data();

    auto order2eval = [&](auto &entry, auto &neighborEntry) {
        const auto &pots = potentials.potentialsOrder2(entry.type);
        auto itPot = pots.find(neighborEntry.type);
        if (itPot != std::end(pots)) {
            Vec3 forceVec{0, 0, 0};
            auto x_ij = bcs::shortestDifference(entry.position(), neighborEntry.position(), box, pbc);
            for (const auto &potential : itPot->second) {
                potential->calculateForceAndEnergy(forceVec, stateModel.energy(), x_ij);
            }
            entry.force += forceVec;
            neighborEntry.force -= forceVec;
            detail::computeVirial<COMPUTE_VIRIAL>(x_ij, forceVec, stateModel.virial());
        }
    };

    auto noop = [](readdy::model::top::GraphTopology* topology){};
    std::vector<readdy::model::top::GraphTopology*> empty;
    readdy::algo::evaluateOnContainers(data, order1eval, neighborList, order2eval, empty, noop);
}

template void MPICalculateForces::performImpl<true>();
template void MPICalculateForces::performImpl<false>();

}