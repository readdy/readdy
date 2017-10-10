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
 * @file SingleCPUCalculateForces.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 11.07.16
 */

#include <readdy/kernel/singlecpu/actions/SCPUCalculateForces.h>

namespace readdy {
namespace kernel {
namespace scpu {
namespace actions {
SCPUCalculateForces::SCPUCalculateForces(SCPUKernel *kernel) : kernel(kernel) {}

void SCPUCalculateForces::perform(const util::PerformanceNode &node) {
    auto t = node.timeit();

    const auto &context = kernel->context();

    auto &stateModel = kernel->getSCPUKernelStateModel();
    auto &data = *stateModel.getParticleData();
    auto &neighborList = *stateModel.getNeighborList();

    stateModel.energy() = 0;
    // update forces and energy order 1 potentials
    {
        const Vec3 zero{0, 0, 0};
        for (auto &e : data) {
            e.force = zero;
            for (const auto &po1 : context.potentials().potentialsOf(e.type)) {
                po1->calculateForceAndEnergy(e.force, stateModel.energy(), e.position());
            }
        }
    }

    // update forces and energy order 2 potentials
    if(!context.potentials().potentialsOrder2().empty()) {
        const auto &difference = context.shortestDifferenceFun();
        Vec3 forceVec{0, 0, 0};
        for (auto it = neighborList.begin(); it != neighborList.end(); ++it) {
            auto i = it->idx1;
            auto j = it->idx2;
            auto &entry_i = data.entry_at(i);
            auto &entry_j = data.entry_at(j);
            const auto &potentials = context.potentials().potentialsOf(entry_i.type, entry_j.type);
            for (const auto &potential : potentials) {
                potential->calculateForceAndEnergy(forceVec, stateModel.energy(), difference(entry_i.position(), entry_j.position()));
                entry_i.force += forceVec;
                entry_j.force += -1 * forceVec;
            }
        }
    }
    // update forces and energy for topologies
    {
        auto taf = kernel->getTopologyActionFactory();
        for ( auto &topology : stateModel.topologies()) {
            if(!topology->isDeactivated()) {
                // calculate bonded potentials
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
            }
        }
    }
}

}
}
}
}