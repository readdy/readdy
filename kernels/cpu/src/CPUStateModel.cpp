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
 * @file CPUStateModel.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 13.07.16
 */

#include <future>
#include <readdy/kernel/cpu/CPUStateModel.h>
#include <readdy/common/thread/barrier.h>
#include <readdy/common/thread/scoped_async.h>
#include <readdy/kernel/cpu/util/config.h>

namespace readdy {
namespace kernel {
namespace cpu {

namespace thd = readdy::util::thread;

using entries_it = CPUStateModel::data_t::entries_t::iterator;
using topologies_it = std::vector<std::unique_ptr<readdy::model::top::GraphTopology>>::const_iterator;
using neighbors_it = decltype(std::declval<readdy::kernel::cpu::model::CPUNeighborList>().cbegin());
using pot1Map = decltype(std::declval<readdy::model::KernelContext>().potentials().potentials_order1());
using pot2Map = decltype(std::declval<readdy::model::KernelContext>().potentials().potentials_order2());
using dist_fun = readdy::model::KernelContext::shortest_dist_fun;
using top_action_factory = readdy::model::top::TopologyActionFactory;

void calculateForcesThread(entries_it begin, entries_it end, neighbors_it neighbors_it,
                           std::promise<double> energyPromise, const CPUStateModel::data_t &data,
                           pot1Map pot1, pot2Map pot2, dist_fun d,
                           topologies_it tbegin, topologies_it tend, top_action_factory const *const taf,
                           const thd::barrier &barrier) {
    double energyUpdate = 0.0;
    for (auto it = begin; it != end; ++it) {
        if (!it->is_deactivated()) {
            readdy::model::Vec3 force{0, 0, 0};
            const auto &myPos = it->position();

            //
            // 1st order potentials
            //
            auto find_it = pot1.find(it->type);
            if (find_it != pot1.end()) {
                for (const auto &potential : find_it->second) {
                    potential->calculateForceAndEnergy(force, energyUpdate, myPos);
                }
            }

            //
            // 2nd order potentials
            //
            double mySecondOrderEnergy = 0.;
            for (const auto neighbor : *neighbors_it) {
                auto &neighborEntry = data.entry_at(neighbor);
                auto potit = pot2.find(std::tie(it->type, neighborEntry.type));
                if (potit != pot2.end()) {
                    auto x_ij = d(myPos, neighborEntry.position());
                    auto distSquared = x_ij * x_ij;
                    for (const auto &potential : potit->second) {
                        if (distSquared < potential->getCutoffRadiusSquared()) {
                            readdy::model::Vec3 updateVec{0, 0, 0};
                            potential->calculateForceAndEnergy(updateVec, mySecondOrderEnergy, x_ij);
                            force += updateVec;
                        }
                    }
                }
            }
            // The contribution of second order potentials must be halved since we parallelise over particles.
            // Thus every particle pair potential is seen twice
            energyUpdate += 0.5 * mySecondOrderEnergy;

            it->force = force;
        }
        ++neighbors_it;
    }

    /*barrier.wait();

    for (auto t_it = tbegin; t_it != tend; ++t_it) {
        const auto &top = *t_it;
        for (const auto &bondedPot : top->getBondedPotentials()) {
            auto energy = bondedPot->createForceAndEnergyAction(taf)->perform();
            energyUpdate += energy;
        }
        for (const auto &anglePot : top->getAnglePotentials()) {
            auto energy = anglePot->createForceAndEnergyAction(taf)->perform();
            energyUpdate += energy;
        }
        for (const auto &torsionPot : top->getTorsionPotentials()) {
            auto energy = torsionPot->createForceAndEnergyAction(taf)->perform();
            energyUpdate += energy;
        }
    }*/

    energyPromise.set_value(energyUpdate);
}

struct CPUStateModel::Impl {
    using particle_t = readdy::model::Particle;
    using reaction_counts_order1_map = std::unordered_map<particle_t::type_type, std::vector<std::size_t>>;
    using reaction_counts_order2_map = std::unordered_map<util::particle_type_pair, std::vector<std::size_t>,
            util::particle_type_pair_hasher, util::particle_type_pair_equal_to>;
    readdy::model::KernelContext *context;
    std::unique_ptr<readdy::kernel::cpu::model::CPUNeighborList> neighborList;
    double currentEnergy = 0;
    std::unique_ptr<readdy::signals::scoped_connection> reorderConnection;
    std::vector<std::unique_ptr<readdy::model::top::GraphTopology>> topologies{};
    top_action_factory const *const topologyActionFactory;
    std::vector<readdy::model::reactions::ReactionRecord> reactionRecords{};
    reaction_counts_order1_map reactionCountsOrder1 {};
    reaction_counts_order2_map reactionCountsOrder2 {};

    template<bool fixpos = true>
    const model::CPUParticleData &cdata() const {
        if (fixpos) particleData->setFixPosFun(context->getFixPositionFun());
        return *particleData;
    }

    template<bool fixpos = true>
    model::CPUParticleData &data() {
        if (fixpos) particleData->setFixPosFun(context->getFixPositionFun());
        return *particleData;
    }

    Impl(readdy::model::KernelContext *context, top_action_factory const *const taf,
         readdy::util::thread::Config const *const config)
            : particleData(std::make_unique<CPUStateModel::data_t>(context)), topologyActionFactory(taf) {
        Impl::context = context;
    }

private:
    std::unique_ptr<readdy::kernel::cpu::model::CPUParticleData> particleData;
};

void CPUStateModel::calculateForces() {
    pimpl->currentEnergy = 0;
    const auto &particleData = pimpl->cdata<true>();
    const auto potOrder1 = pimpl->context->potentials().potentials_order1();
    const auto potOrder2 = pimpl->context->potentials().potentials_order2();
    auto d = pimpl->context->getShortestDifferenceFun();
    {
        std::vector<std::future<double>> energyFutures;
        energyFutures.reserve(config->nThreads());
        {
            std::vector<threading_model> threads;
            threads.reserve(config->nThreads());
            const std::size_t grainSize = (pimpl->cdata().size()) / config->nThreads();
            const std::size_t grainSizeTopologies = pimpl->topologies.size() / config->nThreads();
            auto it_data_end = pimpl->data<false>().end();
            auto it_data = pimpl->data<false>().begin();
            auto it_nl = pimpl->neighborList->begin();
            auto it_tops = pimpl->topologies.cbegin();
            const thd::barrier barrier{config->nThreads()};
            for (auto i = 0; i < config->nThreads() - 1; ++i) {
                std::promise<double> energyPromise;
                energyFutures.push_back(energyPromise.get_future());
                threads.emplace_back(calculateForcesThread, it_data, it_data + grainSize, it_nl,
                                     std::move(energyPromise), std::cref(particleData), potOrder1, potOrder2, d,
                                     it_tops,
                                     it_tops + grainSizeTopologies, pimpl->topologyActionFactory, std::cref(barrier));
                it_nl += grainSize;
                it_data += grainSize;
                it_tops += grainSizeTopologies;
            }
            {
                std::promise<double> lastPromise;
                energyFutures.push_back(lastPromise.get_future());
                threads.emplace_back(calculateForcesThread, it_data, it_data_end, it_nl, std::move(lastPromise),
                                     std::cref(particleData), potOrder1, potOrder2, d, it_tops,
                                     pimpl->topologies.cend(),
                                     pimpl->topologyActionFactory, std::cref(barrier));
            }

        }
        for (auto &f : energyFutures) {
            pimpl->currentEnergy += f.get();
        }
    }

    for (auto t_it = pimpl->topologies.cbegin(); t_it != pimpl->topologies.cend(); ++t_it) {
        const auto &top = *t_it;
        for (const auto &bondedPot : top->getBondedPotentials()) {
            auto action = bondedPot->createForceAndEnergyAction(pimpl->topologyActionFactory);
            auto energy = action->perform();
            pimpl->currentEnergy += energy;
        }
        for (const auto &anglePot : top->getAnglePotentials()) {
            auto energy = anglePot->createForceAndEnergyAction(pimpl->topologyActionFactory)->perform();
            pimpl->currentEnergy += energy;
        }
        for (const auto &torsionPot : top->getTorsionPotentials()) {
            auto energy = torsionPot->createForceAndEnergyAction(pimpl->topologyActionFactory)->perform();
            pimpl->currentEnergy += energy;
        }
    }
}

const std::vector<readdy::model::Vec3> CPUStateModel::getParticlePositions() const {
    const auto &data = pimpl->cdata();
    std::vector<readdy::model::Vec3> target{};
    target.reserve(data.size());
    for (const auto &entry : data) {
        if (!entry.is_deactivated()) target.push_back(entry.position());
    }
    return target;
}

const std::vector<readdy::model::Particle> CPUStateModel::getParticles() const {
    const auto &data = pimpl->cdata();
    std::vector<readdy::model::Particle> result;
    result.reserve(data.size());
    for (const auto &entry : data) {
        if (!entry.is_deactivated()) {
            result.push_back(data.toParticle(entry));
        }
    }
    return result;
}

void CPUStateModel::updateNeighborList() {
    pimpl->neighborList->create();
}

void CPUStateModel::addParticle(const readdy::model::Particle &p) {
    pimpl->data<false>().addParticle(p);
}

void CPUStateModel::addParticles(const std::vector<readdy::model::Particle> &p) {
    pimpl->data<false>().addParticles(p);
}

void CPUStateModel::removeParticle(const readdy::model::Particle &p) {
    pimpl->data<false>().removeParticle(p);
}

double CPUStateModel::getEnergy() const {
    return pimpl->currentEnergy;
}

CPUStateModel::CPUStateModel(readdy::model::KernelContext *const context,
                             readdy::util::thread::Config const *const config,
                             readdy::model::top::TopologyActionFactory const *const taf)
        : pimpl(std::make_unique<Impl>(context, taf, config)), config(config) {
    pimpl->neighborList = std::make_unique<model::CPUNeighborList>(context, *getParticleData(), config);
    pimpl->reorderConnection = std::make_unique<readdy::signals::scoped_connection>(
            pimpl->neighborList->registerReorderEventListener([this](const std::vector<std::size_t> &indices) -> void {
                for (auto &top : pimpl->topologies) {
                    top->permuteIndices(indices);
                }
            }));

}

CPUStateModel::data_t const *const CPUStateModel::getParticleData() const {
    return &pimpl->data<false>();
}

CPUStateModel::data_t *const CPUStateModel::getParticleData() {
    return &pimpl->data<false>();
}

model::CPUNeighborList const *const CPUStateModel::getNeighborList() const {
    return pimpl->neighborList.get();
}

void CPUStateModel::clearNeighborList() {
    pimpl->neighborList->clear();
}

void CPUStateModel::removeAllParticles() {
    pimpl->data<false>().clear();
}

model::CPUNeighborList *const CPUStateModel::getNeighborList() {
    return pimpl->neighborList.get();
}

readdy::model::top::GraphTopology *const
CPUStateModel::addTopology(const std::vector<readdy::model::TopologyParticle> &particles) {
    std::vector<std::size_t> ids = pimpl->data<false>().addTopologyParticles(particles);
    std::vector<particle_type_type> types;
    types.reserve(ids.size());
    for (const auto &p : particles) {
        types.push_back(p.getType());
    }
    pimpl->topologies.push_back(std::make_unique<readdy::model::top::GraphTopology>(std::move(ids), std::move(types),
                                                                                    &pimpl->context->topology_potentials()));
    return pimpl->topologies.back().get();
}

std::vector<readdy::model::reactions::ReactionRecord> &CPUStateModel::reactionRecords() {
    return pimpl->reactionRecords;
}

const std::vector<readdy::model::reactions::ReactionRecord> &CPUStateModel::reactionRecords() const {
    return pimpl->reactionRecords;
}

CPUStateModel::reaction_counts_order1_map &CPUStateModel::reactionCountsOrder1() {
    return pimpl->reactionCountsOrder1;
}

const CPUStateModel::reaction_counts_order1_map &CPUStateModel::reactionCountsOrder1() const {
    return pimpl->reactionCountsOrder1;
}

CPUStateModel::reaction_counts_order2_map &CPUStateModel::reactionCountsOrder2() {
    return pimpl->reactionCountsOrder2;
}

const CPUStateModel::reaction_counts_order2_map &CPUStateModel::reactionCountsOrder2() const {
    return pimpl->reactionCountsOrder2;
}

readdy::model::Particle CPUStateModel::getParticleForIndex(const std::size_t index) const {
    return pimpl->cdata<false>().getParticle(index);
}

void CPUStateModel::expected_n_particles(const std::size_t n) {
    auto& data = pimpl->data<false>();
    if(data.size() < n) {
        data.reserve(n);
    }
}

CPUStateModel::~CPUStateModel() = default;


}
}
}