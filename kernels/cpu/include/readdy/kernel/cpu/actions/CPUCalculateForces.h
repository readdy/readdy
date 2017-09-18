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
    using nl_bounds = std::tuple<nl::NeighborList::const_iterator, nl::NeighborList::const_iterator>;
    using top_bounds = std::tuple<CPUStateModel::topologies_vec::const_iterator, CPUStateModel::topologies_vec::const_iterator>;
public:

    explicit CPUCalculateForces(CPUKernel* kernel) : kernel(kernel) {}

    void perform(const util::PerformanceNode &node) override {
        auto t = node.timeit();

        const auto &context = kernel->getKernelContext();
        const auto &config = kernel->threadConfig();

        auto &stateModel = kernel->getCPUKernelStateModel();
        auto neighborList = stateModel.getNeighborList();
        auto data = neighborList->data();
        auto taf = kernel->getTopologyActionFactory();
        auto &topologies = stateModel.topologies();

        stateModel.energy() = 0;

        const auto &potOrder1 = context.potentials().potentials_order1();
        const auto &potOrder2 = context.potentials().potentials_order2();
        if(!potOrder1.empty() || !potOrder2.empty()) {
            const auto &d = context.shortestDifferenceFun();
            {
                //std::vector<std::future<scalar>> energyFutures;
                //energyFutures.reserve(config->nThreads());
                std::vector<std::promise<scalar>> promises(config.nThreads());
                {
                    const std::size_t grainSize = data->size() / config.nThreads();
                    const std::size_t grainSizeTopologies = topologies.size() / config.nThreads();
                    auto it_data = data->begin();
                    auto it_data_end = data->end();
                    auto it_nl = neighborList->begin();
                    auto it_nl_end = neighborList->end();
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
                        auto nlBounds = std::make_tuple(it_nl, it_nl + grainSize);
                        auto topBounds = std::make_tuple(it_tops, it_tops + grainSizeTopologies);
                        executables.push_back(executor.pack(calculate, dataBounds, nlBounds, topBounds,
                                                            std::ref(promises.at(i)), data, taf, std::cref(context),
                                                            std::cref(barrier)));
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
                                                            std::ref(lastPromise), data, taf, std::cref(context),
                                                            std::cref(barrier)));
                    }
                    executor.execute_and_wait(std::move(executables));

                }
                for (auto &f : promises) {
                    stateModel.energy() += f.get_future().get();
                }
            }
        }

        /*for (auto t_it = topologies.cbegin(); t_it != topologies.cend(); ++t_it) {
            const auto &top = **t_it;
            if (!top.isDeactivated()) {
                for (const auto &bondedPot : top.getBondedPotentials()) {
                    auto action = bondedPot->createForceAndEnergyAction(kernel->getTopologyActionFactory());
                    auto energy = action->perform(&top);
                    stateModel.energy() += energy;
                }
                for (const auto &anglePot : top.getAnglePotentials()) {
                    auto energy = anglePot->createForceAndEnergyAction(kernel->getTopologyActionFactory())->perform(&top);
                    stateModel.energy() += energy;
                }
                for (const auto &torsionPot : top.getTorsionPotentials()) {
                    auto energy = torsionPot->createForceAndEnergyAction(kernel->getTopologyActionFactory())->perform(&top);
                    stateModel.energy() += energy;
                }
            }
        }*/
    }

protected:

    static void calculate(std::size_t /*tid*/, data_bounds dataBounds, nl_bounds nlBounds, top_bounds topBounds,
                   std::promise<scalar>& energyPromise, CPUStateModel::data_type* data,
                   model::top::TopologyActionFactory *taf, const model::KernelContext &context,
                   const readdy::util::thread::barrier &barrier) {
        scalar energyUpdate = 0.0;

        const auto &pot1 = context.potentials().potentials_order1();
        const auto &pot2 = context.potentials().potentials_order2();
        const auto &d = context.shortestDifferenceFun();

        {
            //
            // 1st order potentials
            //
            for (auto it = std::get<0>(dataBounds); it != std::get<1>(dataBounds); ++it) {
                auto &entry = *it;
                if (!entry.deactivated) {
                    Vec3 force{0, 0, 0};
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
            for(auto it = std::get<0>(nlBounds); it != std::get<1>(nlBounds); ++it) {
                auto &entry = data->entry_at(it->current_particle());
                if (!entry.deactivated) {
                    Vec3 force{0, 0, 0};
                    const auto &myPos = entry.pos;


                    auto find_it = pot1.find(entry.type);
                    if (find_it != pot1.end()) {
                        for (const auto &potential : find_it->second) {
                            potential->calculateForceAndEnergy(force, energyUpdate, myPos);
                        }
                    }

                    //
                    // 2nd order potentials
                    //
                    scalar mySecondOrderEnergy = 0.;
                    /*if(!pot2.empty()) */{
                        for(const auto nidx : *it) {
                            auto &neighborEntry = data->entry_at(nidx);
                            auto potit = pot2.find(std::tie(entry.type, neighborEntry.type));
                            if (potit != pot2.end()) {
                                auto x_ij = d(myPos, neighborEntry.pos);
                                auto distSquared = x_ij * x_ij;
                                for (const auto &potential : potit->second) {
                                    if (distSquared < potential->getCutoffRadiusSquared()) {
                                        Vec3 updateVec{0, 0, 0};
                                        potential->calculateForceAndEnergy(updateVec, mySecondOrderEnergy, x_ij);
                                        force += updateVec;
                                    }
                                }
                            }
                        }
                    }
                    // The contribution of second order potentials must be halved since we parallelise over particles.
                    // Thus every particle pair potential is seen twice
                    energyUpdate += 0.5 * mySecondOrderEnergy;

                    entry.force = force;
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

    CPUKernel *const kernel;
};
}
}
}
}
