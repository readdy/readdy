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
    using nl_bounds = std::tuple<nl::CompactCellLinkedList::HEAD::const_iterator, nl::CompactCellLinkedList::HEAD::const_iterator>;
    using top_bounds = std::tuple<CPUStateModel::topologies_vec::const_iterator, CPUStateModel::topologies_vec::const_iterator>;
public:

    explicit CPUCalculateForces(CPUKernel* kernel) : kernel(kernel) {}

    void perform(const util::PerformanceNode &node) override {
        auto t = node.timeit();

        const auto &context = kernel->context();
        const auto &config = kernel->threadConfig();

        auto &stateModel = kernel->getCPUKernelStateModel();
        auto neighborList = stateModel.getNeighborList();
        auto data = stateModel.getParticleData();
        auto taf = kernel->getTopologyActionFactory();
        auto &topologies = stateModel.topologies();

        stateModel.energy() = 0;

        const auto &potOrder1 = context.potentials().potentialsOrder1();
        const auto &potOrder2 = context.potentials().potentialsOrder2();
        if(!potOrder1.empty() || !potOrder2.empty() || !stateModel.topologies().empty()) {
            const auto &d = context.shortestDifferenceFun();
            {
                //std::vector<std::future<scalar>> energyFutures;
                //energyFutures.reserve(config->nThreads());
                std::vector<std::promise<scalar>> promises(config.nThreads());
                {
                    const std::size_t grainSize = data->size() / config.nThreads();
                    const std::size_t grainSizeNeighborList = neighborList->head().size() / config.nThreads();
                    const std::size_t grainSizeTopologies = topologies.size() / config.nThreads();
                    auto it_data = data->begin();
                    auto it_data_end = data->end();
                    auto it_nl = neighborList->head().begin();
                    auto it_nl_end = neighborList->head().end();
                    auto it_tops = topologies.cbegin();
                    auto it_tops_end = topologies.cend();
                    const readdy::util::thread::barrier barrier{config.nThreads()};

                    const auto &executor = *config.executor();
                    std::vector<std::function<void(std::size_t)>> executables;
                    executables.reserve(config.nThreads());

                    for (std::size_t i = 0; i < config.nThreads() - 1; ++i) {
                        //energyFutures.push_back(promises[i].get_future());
                        //ForcesThreadArgs args (, std::cref(barrier));
                        auto dataBounds = std::make_tuple(it_data, it_data + grainSize);
                        auto nlBounds = std::make_tuple(it_nl, it_nl + grainSizeNeighborList);
                        auto topBounds = std::make_tuple(it_tops, it_tops + grainSizeTopologies);
                        executables.push_back(executor.pack(calculate, dataBounds, nlBounds, topBounds,
                                                            std::ref(promises.at(i)), data, std::cref(*neighborList),
                                                            taf, std::cref(context), std::cref(barrier)));
                        it_data = std::get<1>(dataBounds);
                        it_nl = std::get<1>(nlBounds);
                        it_tops = std::get<1>(topBounds);
                    }
                    {
                        std::promise<scalar> &lastPromise = promises.back();
                        //energyFutures.push_back(lastPromise.get_future());
                        //ForcesThreadArgs args (, std::cref(barrier));
                        auto dataBounds = std::make_tuple(it_data, it_data_end);
                        auto nlBounds = std::make_tuple(it_nl, it_nl_end);
                        auto topBounds = std::make_tuple(it_tops, it_tops_end);
                        executables.push_back(executor.pack(calculate, dataBounds, nlBounds, topBounds,
                                                            std::ref(lastPromise), data, std::cref(*neighborList),
                                                            taf, std::cref(context), std::cref(barrier)));
                    }
                    executor.execute_and_wait(std::move(executables));
                }
                for (auto &f : promises) {
                    stateModel.energy() += f.get_future().get();
                }
            }
        }
    }

protected:

    static void calculate(std::size_t /*tid*/, data_bounds dataBounds, nl_bounds nlBounds, top_bounds topBounds,
                          std::promise<scalar>& energyPromise, CPUStateModel::data_type* data,
                          const CPUStateModel::neighbor_list &nl,
                          model::top::TopologyActionFactory *taf, const model::Context &context,
                          const readdy::util::thread::barrier &barrier) {
        scalar energyUpdate = 0.0;

        const auto &pot1 = context.potentials().potentialsOrder1();
        const auto &pot2 = context.potentials().potentialsOrder2();
        const auto &d = context.shortestDifferenceFun();

        {
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
        }

        barrier.wait();

        {
            //
            // 2nd order potentials
            //
            if(!pot2.empty()) {
                for(auto cell = std::get<0>(nlBounds); cell != std::get<1>(nlBounds); ++cell) {
                    for(auto pairIt = nl.cellNeighborsEnd(**cell); pairIt != nl.cellNeighborsEnd(**cell); ++pairIt) {
                        auto pair = *pairIt;
                        auto &entry = data->entry_at(std::get<0>(pair));
                        auto &neighbor = data->entry_at(std::get<1>(pair));
                        if (!entry.deactivated && !neighbor.deactivated) {
                            auto &force = entry.force;
                            const auto &myPos = entry.pos;

                            //
                            // 2nd order potentials
                            //
                            scalar mySecondOrderEnergy = 0.;
                            if(!pot2.empty()) {
                                auto potit = pot2.find(std::tie(entry.type, neighbor.type));
                                if (potit != pot2.end()) {
                                    auto x_ij = d(myPos, neighbor.pos);
                                    auto distSquared = x_ij * x_ij;
                                    for (const auto &potential : potit->second) {
                                        if (distSquared < potential->getCutoffRadiusSquared()) {
                                            Vec3 forceUpdate {0, 0, 0};
                                            potential->calculateForceAndEnergy(forceUpdate, mySecondOrderEnergy, x_ij);
                                            force += forceUpdate;
                                            // log::warn("applying 2nd order potential to particles {} and {} of types {} and {} with resulting force {}", it->current_particle(), nidx, entry.type, neighborEntry.type, forceUpdate);
                                        }
                                    }
                                }
                            }
                            // The contribution of second order potentials must be halved since we parallelize over particles.
                            // Thus every particle pair potential is seen twice
                            energyUpdate += 0.5 * mySecondOrderEnergy;
                        }
                    }

                }
            }

        }

        barrier.wait();

        {
            //
            // external potentials
            //

            for(auto it = std::get<0>(topBounds); it != std::get<1>(topBounds); ++it) {
                const auto &top = *it;
                if(!top->isDeactivated()) {
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
        }

        energyPromise.set_value(energyUpdate);
    }

    CPUKernel *const kernel;
};
}
}
}
}
