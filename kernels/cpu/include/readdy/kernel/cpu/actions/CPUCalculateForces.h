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
 * @file CalculateForces.h
 * @brief << brief description >>
 * @author clonker
 * @date 14.07.16
 */

#pragma once

#include <readdy/model/actions/Actions.h>
#include <readdy/kernel/cpu/CPUKernel.h>
#include <readdy/common/thread/barrier.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace actions {
class CPUCalculateForces : public readdy::model::actions::CalculateForces {
    using data_bounds = std::tuple<data::EntryDataContainer::iterator, data::EntryDataContainer::iterator>;
    using nl_bounds = std::tuple<std::size_t, std::size_t>;
    using top_bounds = std::tuple<CPUStateModel::topologies_vec::const_iterator, CPUStateModel::topologies_vec::const_iterator>;
public:

    explicit CPUCalculateForces(CPUKernel *kernel) : kernel(kernel) {}

    void perform(const util::PerformanceNode &node) override {
        auto t = node.timeit();

        const auto &ctx = kernel->context();

        auto &stateModel = kernel->getCPUKernelStateModel();
        auto neighborList = stateModel.getNeighborList();
        auto data = stateModel.getParticleData();
        auto taf = kernel->getTopologyActionFactory();
        auto &topologies = stateModel.topologies();

        stateModel.energy() = 0;

        const auto &potOrder1 = ctx.potentials().potentialsOrder1();
        const auto &potOrder2 = ctx.potentials().potentialsOrder2();
        if (!potOrder1.empty() || !potOrder2.empty() || !stateModel.topologies().empty()) {
            {
                // todo maybe optimize this by transposing data structure
                auto tClear = node.subnode("clear forces").timeit();
                std::for_each(data->begin(), data->end(), [](auto &entry) {
                    entry.force = {0, 0, 0};
                });
            }
            {
                auto &pool = data->pool();
                std::vector<std::promise<scalar>> promises;
                // 1st order pot + topologies = 2*pool size
                // 2nd order pot <= nl.nCells
                promises.reserve(2*pool.size() + neighborList->nCells());
                {
                    if(!potOrder1.empty()) {
                        // 1st order pot
                        const std::size_t grainSize = data->size() / pool.size();
                        auto it = data->begin();
                        for (auto i = 0_z; i < pool.size()-1; ++i) {
                            auto itNext = std::min(it + grainSize, data->end());
                            if(it != itNext) {
                                promises.emplace_back();
                                auto dataBounds = std::make_tuple(it, itNext);
                                pool.push(calculate_order1, dataBounds, std::ref(promises.back()), data, std::cref(ctx));
                            }
                            it = itNext;
                        }
                        if(it != data->end()) {
                            promises.emplace_back();
                            auto dataBounds = std::make_tuple(it, data->end());
                            pool.push(calculate_order1, dataBounds, std::ref(promises.back()), data, std::cref(ctx));
                        }
                    }
                    if(!topologies.empty()) {
                        const std::size_t grainSize = topologies.size() / pool.size();
                        auto it = topologies.cbegin();
                        for (auto i = 0_z; i < pool.size()-1; ++i) {
                            auto itNext = std::min(it + grainSize, topologies.cend());
                            if(it != itNext) {
                                promises.emplace_back();
                                auto bounds = std::make_tuple(it, itNext);
                                pool.push(calculate_topologies, bounds, taf, std::ref(promises.back()));
                            }
                            it = itNext;
                        }
                        if(it != topologies.cend()) {
                            promises.emplace_back();
                            auto bounds = std::make_tuple(it, topologies.cend());
                            pool.push(calculate_topologies, bounds, taf, std::ref(promises.back()));
                        }
                    }
                    if (!potOrder2.empty()) {
                        for (auto i = 0_z; i < neighborList->nCells(); ++i) {
                            if (!neighborList->cellEmpty(i)) {
                                promises.emplace_back();
                                pool.push(calculate_order2, std::make_tuple(i, i + 1), data, std::cref(*neighborList),
                                          std::ref(promises.back()), std::cref(ctx));
                            }
                        }
                    }
                }
                for (auto &f : promises) {
                    stateModel.energy() += f.get_future().get();
                }
            }
        }
    }

protected:

    static void calculate_order2(std::size_t, nl_bounds nlBounds, CPUStateModel::data_type *data,
                                 const CPUStateModel::neighbor_list &nl, std::promise<scalar> &energyPromise,
                                 const model::Context &context) {
        scalar energyUpdate = 0.0;

        const auto pot2 = context.potentials().potentialsOrder2();
        const auto d = context.shortestDifferenceFun();
        //
        // 2nd order potentials
        //
        for (auto cell = std::get<0>(nlBounds); cell < std::get<1>(nlBounds); ++cell) {
            for (auto particleIt = nl.particlesBegin(cell); particleIt != nl.particlesEnd(cell); ++particleIt) {
                auto &entry = data->entry_at(*particleIt);
                if (entry.deactivated) {
                    log::critical("deactivated particle in neighbor list!");
                    continue;
                }

                nl.forEachNeighbor(*particleIt, cell, [&](auto neighborIndex) {
                    auto &neighbor = data->entry_at(neighborIndex);
                    if (!neighbor.deactivated) {
                        auto &force = entry.force;
                        const auto &myPos = entry.pos;

                        //
                        // 2nd order potentials
                        //
                        scalar mySecondOrderEnergy = 0.;
                        auto potit = pot2.find(std::tie(entry.type, neighbor.type));
                        if (potit != pot2.end()) {
                            auto x_ij = d(myPos, neighbor.pos);
                            auto distSquared = x_ij * x_ij;
                            for (const auto &potential : potit->second) {
                                if (distSquared < potential->getCutoffRadiusSquared()) {
                                    Vec3 forceUpdate{0, 0, 0};
                                    potential->calculateForceAndEnergy(forceUpdate, mySecondOrderEnergy, x_ij);
                                    force += forceUpdate;
                                }
                            }
                        }
                        // The contribution of second order potentials must be halved since we parallelize over particles.
                        // Thus every particle pair potential is seen twice
                        energyUpdate += 0.5 * mySecondOrderEnergy;
                    } else {
                        log::critical("disabled neighbour");
                    }
                });
            }

        }
        energyPromise.set_value(energyUpdate);

    }

    static void calculate_topologies(std::size_t /*tid*/, top_bounds topBounds, model::top::TopologyActionFactory *taf,
                                     std::promise<scalar> &energyPromise) {
        scalar energyUpdate = 0.0;
        for (auto it = std::get<0>(topBounds); it != std::get<1>(topBounds); ++it) {
            const auto &top = *it;
            if (!top->isDeactivated()) {
                for (const auto &bondedPot : top->getBondedPotentials()) {
                    auto energy = bondedPot->createForceAndEnergyAction(taf)->perform(top.get());
                    energyUpdate += energy;
                }
                for (const auto &anglePot : top->getAnglePotentials()) {
                    auto energy = anglePot->createForceAndEnergyAction(taf)->perform(top.get());
                    energyUpdate += energy;
                }
                for (const auto &torsionPot : top->getTorsionPotentials()) {
                    auto energy = torsionPot->createForceAndEnergyAction(taf)->perform(top.get());
                    energyUpdate += energy;
                }
            }
        }

        energyPromise.set_value(energyUpdate);
    }


    static void calculate_order1(std::size_t /*tid*/, data_bounds dataBounds,
                                 std::promise<scalar> &energyPromise, CPUStateModel::data_type *data,
                                 const model::Context &context) {
        scalar energyUpdate = 0.0;

        auto pot1 = context.potentials().potentialsOrder1();
        const auto d = context.shortestDifferenceFun();

        //
        // 1st order potentials
        //
        for (auto it = std::get<0>(dataBounds); it != std::get<1>(dataBounds); ++it) {
            auto &entry = *it;
            if (!entry.deactivated) {
                auto &force = entry.force;
                force = {c_::zero, c_::zero, c_::zero};
                const auto &myPos = entry.pos;
                auto find_it = pot1.find(entry.type);
                if (find_it != pot1.end()) {
                    for (const auto &potential : find_it->second) {
                        potential->calculateForceAndEnergy(force, energyUpdate, myPos);
                    }
                }
            }
        }
        energyPromise.set_value(energyUpdate);
    }

    CPUKernel *const kernel;
};
}
}
}
}
