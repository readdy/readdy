/********************************************************************
 * Copyright © 2018 Computational Molecular Biology Group,          *
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
 * @file CPUCalculateForces.cpp
 * @author clonker
 * @date 1/3/18
 */

#include "readdy/kernel/cpu/actions/CPUCalculateForces.h"

namespace readdy {
namespace kernel {
namespace cpu {
namespace actions {

void CPUCalculateForces::perform(const util::PerformanceNode &node) {
    auto t = node.timeit();

    const auto &ctx = kernel->context();

    auto &stateModel = kernel->getCPUKernelStateModel();
    auto neighborList = stateModel.getNeighborList();
    auto data = stateModel.getParticleData();
    auto taf = kernel->getTopologyActionFactory();
    auto &topologies = stateModel.topologies();

    stateModel.energy() = 0;
    stateModel.virial() = Matrix33{{{0, 0, 0, 0, 0, 0, 0, 0, 0}}};

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
            std::vector<std::promise<Matrix33>> virialPromises;
            // 1st order pot + topologies = 2*pool size
            // 2nd order pot <= nl.nCells
            size_t nThreads = pool.size();
            auto numberTasks = (!potOrder1.empty() ? nThreads : 0)
                               + (!potOrder2.empty() ? nThreads : 0)
                               + (!topologies.empty() ? nThreads : 0);
            {
                const auto &nTasks = node.subnode("create tasks");
                auto tTasks = nTasks.timeit();
                size_t nCells = neighborList->nCells();
                promises.reserve(numberTasks);
                virialPromises.reserve(!potOrder2.empty() ? nThreads : 0);
                if (!potOrder1.empty()) {
                    // 1st order pot
                    auto tO1 = nTasks.subnode("order1").timeit();
                    std::vector<std::function<void(std::size_t)>> tasks;
                    tasks.reserve(nThreads);

                    const std::size_t grainSize = data->size() / nThreads;
                    auto it = data->begin();
                    for (auto i = 0_z; i < nThreads - 1; ++i) {
                        auto itNext = std::min(it + grainSize, data->end());
                        if (it != itNext) {
                            promises.emplace_back();
                            auto dataBounds = std::make_tuple(it, itNext);
                            tasks.push_back(pool.pack(calculate_order1, dataBounds, std::ref(promises.back()), data,
                                                      ctx.potentials().potentialsOrder1()));
                        }
                        it = itNext;
                    }
                    if (it != data->end()) {
                        promises.emplace_back();
                        auto dataBounds = std::make_tuple(it, data->end());
                        tasks.push_back(pool.pack(calculate_order1, dataBounds, std::ref(promises.back()), data,
                                                  ctx.potentials().potentialsOrder1()));
                    }
                    {
                        auto tPush = nTasks.subnode("execute order 1 tasks and wait").timeit();
                        auto futures = pool.pushAll(std::move(tasks));
                        std::vector<util::thread::joining_future<void>> joiningFutures;
                        std::transform(futures.begin(), futures.end(), std::back_inserter(joiningFutures),
                                       [](auto &&future) {
                                           return util::thread::joining_future<void>{std::move(future)};
                                       });
                    }
                }
                if (!topologies.empty()) {
                    auto tTops = nTasks.subnode("topologies").timeit();
                    std::vector<std::function<void(std::size_t)>> tasks;
                    tasks.reserve(nThreads);
                    const std::size_t grainSize = topologies.size() / nThreads;
                    auto it = topologies.cbegin();
                    for (auto i = 0_z; i < nThreads - 1; ++i) {
                        auto itNext = std::min(it + grainSize, topologies.cend());
                        if (it != itNext) {
                            promises.emplace_back();
                            auto bounds = std::make_tuple(it, itNext);
                            tasks.push_back(
                                    pool.pack(calculate_topologies, bounds, taf, std::ref(promises.back())));
                        }
                        it = itNext;
                    }
                    if (it != topologies.cend()) {
                        promises.emplace_back();
                        auto bounds = std::make_tuple(it, topologies.cend());
                        tasks.push_back(pool.pack(calculate_topologies, bounds, taf, std::ref(promises.back())));
                    }
                    {
                        auto tPush = nTasks.subnode("execute topology tasks and wait").timeit();
                        auto futures = pool.pushAll(std::move(tasks));
                        std::vector<util::thread::joining_future<void>> joiningFutures;
                        std::transform(futures.begin(), futures.end(), std::back_inserter(joiningFutures),
                                       [](auto &&future) {
                                           return util::thread::joining_future<void>{std::move(future)};
                                       });
                    }
                }
                if (!potOrder2.empty()) {
                    auto tO2 = nTasks.subnode("order2").timeit();
                    std::vector<std::function<void(std::size_t)>> tasks;
                    tasks.reserve(nThreads);
                    auto granularity = nThreads;
                    const std::size_t grainSize = nCells / granularity;
                    auto it = 0_z;
                    for (auto i = 0_z; i < granularity - 1; ++i) {
                        auto itNext = std::min(it + grainSize, nCells);
                        if (it != itNext) {
                            promises.emplace_back();
                            virialPromises.emplace_back();
                            if(ctx.recordVirial()) {
                                tasks.push_back(pool.pack(
                                        calculate_order2<true>, std::make_tuple(it, itNext), data,
                                        std::cref(*neighborList), std::ref(promises.back()),
                                        std::ref(virialPromises.back()), ctx.potentials().potentialsOrder2(),
                                        ctx.boxSize(), ctx.periodicBoundaryConditions()
                                ));
                            } else {
                                tasks.push_back(pool.pack(
                                        calculate_order2<false>, std::make_tuple(it, itNext), data,
                                        std::cref(*neighborList), std::ref(promises.back()),
                                        std::ref(virialPromises.back()), ctx.potentials().potentialsOrder2(),
                                        ctx.boxSize(), ctx.periodicBoundaryConditions()
                                ));
                            }
                        }
                        it = itNext;
                    }
                    if (it != nCells) {
                        promises.emplace_back();
                        virialPromises.emplace_back();
                        if(ctx.recordVirial()) {
                            tasks.push_back(pool.pack(
                                    calculate_order2<true>, std::make_tuple(it, nCells), data, std::cref(*neighborList),
                                    std::ref(promises.back()), std::ref(virialPromises.back()),
                                    ctx.potentials().potentialsOrder2(), ctx.boxSize(), ctx.periodicBoundaryConditions()
                            ));
                        } else {
                            tasks.push_back(pool.pack(
                                    calculate_order2<false>, std::make_tuple(it, nCells), data, std::cref(*neighborList),
                                    std::ref(promises.back()), std::ref(virialPromises.back()),
                                    ctx.potentials().potentialsOrder2(), ctx.boxSize(), ctx.periodicBoundaryConditions()
                            ));
                        }
                    }
                    {
                        auto tPush = nTasks.subnode("execute order 2 tasks and wait").timeit();
                        auto futures = pool.pushAll(std::move(tasks));
                        std::vector<util::thread::joining_future<void>> joiningFutures;
                        std::transform(futures.begin(), futures.end(), std::back_inserter(joiningFutures),
                                       [](auto &&future) {
                                           return util::thread::joining_future<void>{std::move(future)};
                                       });
                    }
                }
            }

            {
                auto tFutures = node.subnode("get energy futures").timeit();
                for (auto &f : promises) {
                    stateModel.energy() += f.get_future().get();
                }
            }
            {
                auto tVirialFutures = node.subnode("get virial futures").timeit();
                for (auto &f : virialPromises) {
                    stateModel.virial() += f.get_future().get();
                }
            }
        }
    }
}

template<bool COMPUTE_VIRIAL>
void CPUCalculateForces::calculate_order2(std::size_t, nl_bounds nlBounds,
                                          CPUStateModel::data_type *data, const CPUStateModel::neighbor_list &nl,
                                          std::promise<scalar> &energyPromise, std::promise<Matrix33> &virialPromise,
                                          model::potentials::PotentialRegistry::PotentialsO2Map pot2,
                                          model::Context::BoxSize box, model::Context::PeriodicBoundaryConditions pbc) {
    scalar energyUpdate = 0.0;
    Matrix33 virialUpdate{{{0, 0, 0, 0, 0, 0, 0, 0, 0}}};

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
                        auto x_ij = bcs::shortestDifference(myPos, neighbor.pos, box.data(), pbc.data());
                        auto distSquared = x_ij * x_ij;
                        for (const auto &potential : potit->second) {
                            if (distSquared < potential->getCutoffRadiusSquared()) {
                                Vec3 forceUpdate{0, 0, 0};
                                potential->calculateForceAndEnergy(forceUpdate, mySecondOrderEnergy, x_ij);
                                force += forceUpdate;
                                if(COMPUTE_VIRIAL && *particleIt < neighborIndex) {
                                    virialUpdate += math::outerProduct<Matrix33>(-1.*x_ij, forceUpdate);
                                }
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
    virialPromise.set_value(virialUpdate);

}

void CPUCalculateForces::calculate_topologies(std::size_t, top_bounds topBounds,
                                              model::top::TopologyActionFactory *taf,
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

void CPUCalculateForces::calculate_order1(std::size_t, data_bounds dataBounds,
                                          std::promise<scalar> &energyPromise, CPUStateModel::data_type *data,
                                          model::potentials::PotentialRegistry::PotentialsO1Map pot1) {
    scalar energyUpdate = 0.0;

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
}
}
}
}
