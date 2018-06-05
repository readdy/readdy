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
 * @file SCPUCalculateForces.h
 * @brief << brief description >>
 * @author clonker
 * @date 20.06.16
 */
#pragma once
#include <readdy/model/actions/Actions.h>
#include <readdy/kernel/singlecpu/SCPUKernel.h>
#include <readdy/common/boundary_condition_operations.h>
#include <readdy/common/algorithm.h>

namespace readdy {
namespace kernel {
namespace scpu {
namespace actions {

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

class SCPUCalculateForces : public readdy::model::actions::CalculateForces {
public:
    explicit SCPUCalculateForces(SCPUKernel *kernel) : kernel(kernel) {};

    void perform(const util::PerformanceNode &node) override {
        const auto &context = kernel->context();
        if(context.recordVirial()) {
            performImpl<true>(node);
        } else {
            performImpl<false>(node);
        }

    };

private:

    template<bool COMPUTE_VIRIAL>
    void performImpl(const util::PerformanceNode &node) {
        auto t = node.timeit();

        const auto &context = kernel->context();

        auto &stateModel = kernel->getSCPUKernelStateModel();
        auto &data = *stateModel.getParticleData();
        auto &neighborList = *stateModel.getNeighborList();

        stateModel.energy() = 0;
        stateModel.virial() = Matrix33{{{0, 0, 0, 0, 0, 0, 0, 0, 0}}};

        const auto &potentials = context.potentials();
        auto &topologies = stateModel.topologies();
        if (!potentials.potentialsOrder1().empty() || !potentials.potentialsOrder2().empty() || !topologies.empty()) {
            auto tClear = node.subnode("clear").timeit();
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


        // topology eval
        auto taf = kernel->getTopologyActionFactory();
        auto topologyEval = [&](auto &topology){
            for (const auto &bondedPot : topology->getBondedPotentials()) {
                auto energy = bondedPot->createForceAndEnergyAction(taf)->perform(topology.get());
                stateModel.energy() += energy;
            }
            for (const auto &anglePot : topology->getAnglePotentials()) {
                auto energy = anglePot->createForceAndEnergyAction(taf)->perform(topology.get());
                stateModel.energy() += energy;
            }
            for (const auto &torsionPot : topology->getTorsionPotentials()) {
                auto energy = torsionPot->createForceAndEnergyAction(taf)->perform(topology.get());
                stateModel.energy() += energy;
            }
        };

        readdy::algo::evaluateOnContainers(data, order1eval, neighborList, order2eval, topologies, topologyEval, node);
    }
    SCPUKernel *kernel;
};
}
}
}
}
