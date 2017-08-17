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
#include <readdy/kernel/cpu/nl/NeighborList.h>
#include <readdy/common/index_persistent_vector.h>

namespace readdy {
namespace kernel {
namespace cpu {

namespace thd = readdy::util::thread;

using entries_it = CPUStateModel::data_t::entries_t::iterator;
using topologies_it = std::vector<std::unique_ptr<readdy::model::top::GraphTopology>>::const_iterator;
using pot1Map = readdy::model::potentials::PotentialRegistry::potential_o1_registry;
using pot2Map = readdy::model::potentials::PotentialRegistry::potential_o2_registry;
using dist_fun = readdy::model::KernelContext::shortest_dist_fun;
using top_action_factory = readdy::model::top::TopologyActionFactory;

void calculateForcesThread(std::size_t, entries_it begin, entries_it end, neighbor_list::const_iterator neighbors_it,
                           std::promise<scalar>& energyPromise, const CPUStateModel::data_t& data, const pot1Map& pot1,
                           const pot2Map& pot2, const dist_fun& d) {
    scalar energyUpdate = 0.0;
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
            scalar mySecondOrderEnergy = 0.;
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
    // todo reactivate this ?
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
    using reaction_counts_order1_map = CPUStateModel::reaction_counts_order1_map;
    using reaction_counts_order2_map = CPUStateModel::reaction_counts_order2_map;
    readdy::model::KernelContext *context;
    std::unique_ptr<neighbor_list> neighborList;
    bool initial_neighbor_list_setup {true};
    scalar currentEnergy = 0;
    std::unique_ptr<readdy::signals::scoped_connection> reorderConnection;
    topologies_vec topologies{};
    top_action_factory const *const topologyActionFactory;
    std::vector<readdy::model::reactions::ReactionRecord> reactionRecords{};
    std::pair<reaction_counts_order1_map, reaction_counts_order2_map> reactionCounts;

    const model::CPUParticleData &cdata() const {
        return *particleData;
    }

    model::CPUParticleData &data() {
        return *particleData;
    }

    Impl(readdy::model::KernelContext *context, top_action_factory const *const taf,
         readdy::util::thread::Config const *const config)
            : particleData(std::make_unique<CPUStateModel::data_t>(context, *config)), topologyActionFactory(taf) {
        Impl::context = context;
    }

private:
    std::unique_ptr<readdy::kernel::cpu::model::CPUParticleData> particleData;
};

void CPUStateModel::calculateForces() {
    pimpl->currentEnergy = 0;
    const auto &particleData = pimpl->cdata();
    const auto &potOrder1 = pimpl->context->potentials().potentials_order1();
    const auto &potOrder2 = pimpl->context->potentials().potentials_order2();
    const auto &d = pimpl->context->getShortestDifferenceFun();
    {
        //std::vector<std::future<scalar>> energyFutures;
        //energyFutures.reserve(config->nThreads());
        std::vector<std::promise<scalar>> promises (config->nThreads());
        {
            const std::size_t grainSize = (pimpl->cdata().size()) / config->nThreads();
            const std::size_t grainSizeTopologies = pimpl->topologies.size() / config->nThreads();
            auto it_data_end = pimpl->data().end();
            auto it_data = pimpl->data().begin();
            auto it_nl = pimpl->neighborList->begin();
            auto it_tops = pimpl->topologies.cbegin();
            const thd::barrier barrier{config->nThreads()};

            const auto& executor = *config->executor();
            std::vector<std::function<void(std::size_t)>> executables;
            executables.reserve(config->nThreads());

            for (std::size_t i = 0; i < config->nThreads() - 1; ++i) {
                //energyFutures.push_back(promises[i].get_future());
                //ForcesThreadArgs args (, std::cref(barrier));
                executables.push_back(executor.pack(calculateForcesThread, it_data, it_data + grainSize, it_nl, std::ref(promises.at(i)),
                                                    std::cref(particleData), std::cref(potOrder1), std::cref(potOrder2),
                                                    std::cref(d)));
                it_nl += grainSize;
                it_data += grainSize;
                it_tops += grainSizeTopologies;
            }
            {
                std::promise<scalar> &lastPromise = promises.back();
                //energyFutures.push_back(lastPromise.get_future());
                //ForcesThreadArgs args (, std::cref(barrier));
                executables.push_back(
                        executor.pack(calculateForcesThread, it_data, it_data_end, it_nl, std::ref(lastPromise),
                                      std::cref(particleData), std::cref(potOrder1), std::cref(potOrder2),
                                      std::cref(d)));
            }
            executor.execute_and_wait(std::move(executables));

        }
        for (auto &f : promises) {
            pimpl->currentEnergy += f.get_future().get();
        }
    }

    for (auto t_it = pimpl->topologies.cbegin(); t_it != pimpl->topologies.cend(); ++t_it) {
        const auto &top = **t_it;
        if (!top.isDeactivated()) {
            for (const auto &bondedPot : top.getBondedPotentials()) {
                auto action = bondedPot->createForceAndEnergyAction(pimpl->topologyActionFactory);
                auto energy = action->perform(&top);
                pimpl->currentEnergy += energy;
            }
            for (const auto &anglePot : top.getAnglePotentials()) {
                auto energy = anglePot->createForceAndEnergyAction(pimpl->topologyActionFactory)->perform(&top);
                pimpl->currentEnergy += energy;
            }
            for (const auto &torsionPot : top.getTorsionPotentials()) {
                auto energy = torsionPot->createForceAndEnergyAction(pimpl->topologyActionFactory)->perform(&top);
                pimpl->currentEnergy += energy;
            }
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
    if(pimpl->initial_neighbor_list_setup) {
        pimpl->neighborList->set_up();
        pimpl->initial_neighbor_list_setup = false;
    } else {
        pimpl->neighborList->update();
    }
}

void CPUStateModel::addParticle(const readdy::model::Particle &p) {
    pimpl->data().addParticle(p);
}

void CPUStateModel::addParticles(const std::vector<readdy::model::Particle> &p) {
    pimpl->data().addParticles(p);
}

void CPUStateModel::removeParticle(const readdy::model::Particle &p) {
    pimpl->data().removeParticle(p);
}

scalar CPUStateModel::getEnergy() const {
    return pimpl->currentEnergy;
}

CPUStateModel::CPUStateModel(readdy::model::KernelContext *const context,
                             readdy::util::thread::Config const *const config,
                             readdy::model::top::TopologyActionFactory const *const taf)
        : pimpl(std::make_unique<Impl>(context, taf, config)), config(config) {
    pimpl->neighborList = std::make_unique<neighbor_list>(*getParticleData(), *context, *config);
    pimpl->reorderConnection = std::make_unique<readdy::signals::scoped_connection>(
            pimpl->data().registerReorderEventListener([this](const std::vector<std::size_t> &indices) -> void {
                for (auto &top : pimpl->topologies) {
                    if(!top->isDeactivated()) top->permuteIndices(indices);
                }
            }));
}

CPUStateModel::data_t const *const CPUStateModel::getParticleData() const {
    return &pimpl->data();
}

CPUStateModel::data_t *const CPUStateModel::getParticleData() {
    return &pimpl->data();
}

neighbor_list const *const CPUStateModel::getNeighborList() const {
    return pimpl->neighborList.get();
}

void CPUStateModel::clearNeighborList() {
    pimpl->neighborList->clear();
}

void CPUStateModel::removeAllParticles() {
    pimpl->data().clear();
}

neighbor_list *const CPUStateModel::getNeighborList() {
    return pimpl->neighborList.get();
}

readdy::model::top::GraphTopology *const
CPUStateModel::addTopology(const std::vector<readdy::model::TopologyParticle> &particles) {
    std::vector<std::size_t> ids = pimpl->data().addTopologyParticles(particles);
    std::vector<particle_type_type> types;
    types.reserve(ids.size());
    for (const auto &p : particles) {
        types.push_back(p.getType());
    }
    auto it = pimpl->topologies.push_back(
            std::make_unique<topology>(std::move(ids), std::move(types), pimpl->context->topology_potentials())
    );
    const auto idx = std::distance(topologies().begin(), it);
    for(const auto p : (*it)->getParticles()) {
        pimpl->data().entry_at(p).topology_index = idx;
    }
    return it->get();
}

std::vector<readdy::model::reactions::ReactionRecord> &CPUStateModel::reactionRecords() {
    return pimpl->reactionRecords;
}

const std::vector<readdy::model::reactions::ReactionRecord> &CPUStateModel::reactionRecords() const {
    return pimpl->reactionRecords;
}

readdy::model::Particle CPUStateModel::getParticleForIndex(const std::size_t index) const {
    return pimpl->cdata().getParticle(index);
}

void CPUStateModel::expected_n_particles(const std::size_t n) {
    auto& data = pimpl->data();
    if(data.size() < n) {
        data.reserve(n);
    }
}

const std::pair<CPUStateModel::reaction_counts_order1_map, CPUStateModel::reaction_counts_order2_map> &CPUStateModel::reactionCounts() const {
    return pimpl->reactionCounts;
}

std::pair<CPUStateModel::reaction_counts_order1_map, CPUStateModel::reaction_counts_order2_map> &CPUStateModel::reactionCounts() {
    return pimpl->reactionCounts;
}

const CPUStateModel::topologies_vec &CPUStateModel::topologies() const {
    return pimpl->topologies;
}

CPUStateModel::topologies_vec &CPUStateModel::topologies() {
    return pimpl->topologies;
}

std::vector<readdy::model::top::GraphTopology const*> CPUStateModel::getTopologies() const {
    std::vector<readdy::model::top::GraphTopology const*> result;
    result.reserve(pimpl->topologies.size() - pimpl->topologies.n_deactivated());
    for(const auto& top : pimpl->topologies) {
        if(!top->isDeactivated()) {
            result.push_back(top.get());
        }
    }
    return result;
}

particle_type_type CPUStateModel::getParticleType(const std::size_t index) const {
    return pimpl->cdata().entry_at(index).type;
}

const readdy::model::top::GraphTopology *CPUStateModel::getTopologyForParticle(readdy::model::top::Topology::particle_t particle) const {
    const auto& entry = pimpl->cdata().entry_at(particle);
    if(!entry.deactivated) {
        if(entry.topology_index >= 0) {
            return pimpl->topologies.at(static_cast<topologies_vec::size_type>(entry.topology_index)).get();
        }
        log::trace("requested particle {} of type {} had no assigned topology", particle, entry.type);
        return nullptr;
    }
    throw std::logic_error(fmt::format("requested particle was deactivated in getTopologyForParticle(p={})", particle));
}

readdy::model::top::GraphTopology *CPUStateModel::getTopologyForParticle(readdy::model::top::Topology::particle_t particle) {
    const auto& entry = pimpl->data().entry_at(particle);
    if(!entry.deactivated) {
        if(entry.topology_index >= 0) {
            return pimpl->topologies.at(static_cast<topologies_vec::size_type>(entry.topology_index)).get();
        }
        log::trace("requested particle {} of type {} had no assigned topology", particle, entry.type);
        return nullptr;
    }
    throw std::logic_error(fmt::format("requested particle was deactivated in getTopologyForParticle(p={})", particle));
}

CPUStateModel::~CPUStateModel() = default;


}
}
}