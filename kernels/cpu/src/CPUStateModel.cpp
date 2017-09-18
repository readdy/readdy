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
#include <readdy/kernel/cpu/nl/AdaptiveNeighborList.h>
#include <readdy/kernel/cpu/nl/CellDecompositionNeighborList.h>
#include <readdy/kernel/cpu/nl/CLLNeighborList.h>

namespace readdy {
namespace kernel {
namespace cpu {

namespace thd = readdy::util::thread;

using entries_it = CPUStateModel::data_type::Entries::iterator;
using topologies_it = std::vector<std::unique_ptr<readdy::model::top::GraphTopology>>::const_iterator;
using pot1Map = readdy::model::potentials::PotentialRegistry::potential_o1_registry;
using pot2Map = readdy::model::potentials::PotentialRegistry::potential_o2_registry;
using dist_fun = readdy::model::KernelContext::shortest_dist_fun;
using top_action_factory = readdy::model::top::TopologyActionFactory;

void calculateForcesThread(std::size_t /*tid*/, const neighbor_list::const_iterator &begin,
                           const neighbor_list::const_iterator &end,
                           std::promise<scalar>& energyPromise, CPUStateModel::data_type* data, const pot1Map& pot1,
                           const pot2Map& pot2, const dist_fun& d) {
    scalar energyUpdate = 0.0;
    for (auto it = begin; it != end; ++it) {
        auto &entry = data->entry_at(it->current_particle());
        if (!entry.deactivated) {
            Vec3 force{0, 0, 0};
            const auto &myPos = entry.pos;

            //
            // 1st order potentials
            //
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
    scalar currentEnergy = 0;
    std::unique_ptr<readdy::signals::scoped_connection> reorderConnection;
    topologies_vec topologies{};
    top_action_factory const *const topologyActionFactory;
    std::vector<readdy::model::reactions::ReactionRecord> reactionRecords{};
    std::pair<reaction_counts_order1_map, reaction_counts_order2_map> reactionCounts;

    Impl(readdy::model::KernelContext *context, top_action_factory const *const taf,
         readdy::util::thread::Config const *const config) : topologyActionFactory(taf), context(context) {}
};

const std::vector<Vec3> CPUStateModel::getParticlePositions() const {
    const auto &data = *pimpl->neighborList->data();
    std::vector<Vec3> target{};
    target.reserve(data.size());
    for (const auto &entry : data) {
        if (!entry.deactivated) target.push_back(entry.pos);
    }
    return target;
}

const std::vector<readdy::model::Particle> CPUStateModel::getParticles() const {
    const auto &data = *pimpl->neighborList->data();
    std::vector<readdy::model::Particle> result;
    result.reserve(data.size());
    for (const auto &entry : data) {
        if (!entry.deactivated) {
            result.push_back(data.toParticle(entry));
        }
    }
    return result;
}

void CPUStateModel::updateNeighborList() {
    updateNeighborList({});
}

void CPUStateModel::addParticle(const readdy::model::Particle &p) {
    getParticleData()->addParticle(p);
}

void CPUStateModel::addParticles(const std::vector<readdy::model::Particle> &p) {
    getParticleData()->addParticles(p);
}

void CPUStateModel::removeParticle(const readdy::model::Particle &p) {
    getParticleData()->removeParticle(p);
}

CPUStateModel::CPUStateModel(readdy::model::KernelContext *const context,
                             readdy::util::thread::Config const *const config,
                             readdy::model::top::TopologyActionFactory const *const taf)
        : pimpl(std::make_unique<Impl>(context, taf, config)), config(config) {
    pimpl->neighborList = std::unique_ptr<neighbor_list>(
            // new nl::CellDecompositionNeighborList(*getParticleData(), *pimpl->context, *config)
            //new nl::CompactCLLNeighborList(1, *pimpl->context, *config)
            new nl::DynamicCLLNeighborList(*pimpl->context, *config)
    );
    pimpl->reorderConnection = std::make_unique<readdy::signals::scoped_connection>(
            getParticleData()->registerReorderEventListener([this](const std::vector<std::size_t> &indices) -> void {
                for (auto &top : pimpl->topologies) {
                    if(!top->isDeactivated()) top->permuteIndices(indices);
                }
            }));
}

CPUStateModel::data_type const *const CPUStateModel::getParticleData() const {
    return pimpl->neighborList->data();
}

CPUStateModel::data_type *const CPUStateModel::getParticleData() {
    return pimpl->neighborList->data();
}

neighbor_list const *const CPUStateModel::getNeighborList() const {
    return pimpl->neighborList.get();
}

void CPUStateModel::clearNeighborList() {
    clearNeighborList({});
}

void CPUStateModel::removeAllParticles() {
    getParticleData()->clear();
}

neighbor_list *const CPUStateModel::getNeighborList() {
    return pimpl->neighborList.get();
}

readdy::model::top::GraphTopology *const
CPUStateModel::addTopology(topology_type_type type, const std::vector<readdy::model::TopologyParticle> &particles) {
    std::vector<std::size_t> ids = getParticleData()->addTopologyParticles(particles);
    std::vector<particle_type_type> types;
    types.reserve(ids.size());
    for (const auto &p : particles) {
        types.push_back(p.getType());
    }
    auto it = pimpl->topologies.push_back(
            std::make_unique<topology>(type, std::move(ids), std::move(types),
                                       pimpl->context->topology_registry().potential_configuration())
    );
    const auto idx = std::distance(topologies().begin(), it);
    for(const auto p : (*it)->getParticles()) {
        getParticleData()->entry_at(p).topology_index = idx;
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
    return getParticleData()->getParticle(index);
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

std::vector<readdy::model::top::GraphTopology*> CPUStateModel::getTopologies() {
    std::vector<readdy::model::top::GraphTopology*> result;
    result.reserve(pimpl->topologies.size() - pimpl->topologies.n_deactivated());
    for(const auto& top : pimpl->topologies) {
        if(!top->isDeactivated()) {
            result.push_back(top.get());
        }
    }
    return result;
}

particle_type_type CPUStateModel::getParticleType(const std::size_t index) const {
    return getParticleData()->entry_at(index).type;
}

const readdy::model::top::GraphTopology *CPUStateModel::getTopologyForParticle(readdy::model::top::Topology::particle_index particle) const {
    const auto& entry = getParticleData()->entry_at(particle);
    if(!entry.deactivated) {
        if(entry.topology_index >= 0) {
            return pimpl->topologies.at(static_cast<topologies_vec::size_type>(entry.topology_index)).get();
        }
        log::trace("requested particle {} of type {} had no assigned topology", particle, entry.type);
        return nullptr;
    }
    throw std::logic_error(fmt::format("requested particle was deactivated in getTopologyForParticle(p={})", particle));
}

readdy::model::top::GraphTopology *CPUStateModel::getTopologyForParticle(readdy::model::top::Topology::particle_index particle) {
    const auto& entry = getParticleData()->entry_at(particle);
    if(!entry.deactivated) {
        if(entry.topology_index >= 0) {
            return pimpl->topologies.at(static_cast<topologies_vec::size_type>(entry.topology_index)).get();
        }
        log::trace("requested particle {} of type {} had no assigned topology", particle, entry.type);
        return nullptr;
    }
    throw std::logic_error(fmt::format("requested particle was deactivated in getTopologyForParticle(p={})", particle));
}

void CPUStateModel::insert_topology(CPUStateModel::topology &&top) {
    auto& topologies = pimpl->topologies;
    auto it = topologies.push_back(std::make_unique<topology>(std::move(top)));
    auto idx = std::distance(topologies.begin(), it);
    const auto& particles = it->get()->getParticles();
    auto& data = *getParticleData();
    std::for_each(particles.begin(), particles.end(), [idx, &data](const topology::particle_index p) {
        data.entry_at(p).topology_index = idx;
    });
}

void CPUStateModel::updateNeighborList(const util::PerformanceNode &node) {
    pimpl->neighborList->update(node.subnode("update"));
}

void CPUStateModel::clearNeighborList(const util::PerformanceNode &node) {
    pimpl->neighborList->clear(node.subnode("clear"));
}

void CPUStateModel::initializeNeighborList(scalar skin) {
    initializeNeighborList(skin, {});
}

void CPUStateModel::initializeNeighborList(scalar skin, const util::PerformanceNode &node) {
    pimpl->neighborList->skin() = skin;
    pimpl->neighborList->set_up(node.subnode("set_up"));
}

void CPUStateModel::configure(const readdy::conf::cpu::Configuration &configuration) {
    const auto& nl = configuration.neighborList;

    if(nl.type == "CellDecomposition") {
        (*pimpl).neighborList = std::unique_ptr<neighbor_list>(
                new nl::CellDecompositionNeighborList(getParticleData(), *pimpl->context, *config)
        );
        log::debug("Using CPU/CellDecomposition neighborList");
    } else if(nl.type == "ContiguousCLL") {
        (*pimpl).neighborList = std::unique_ptr<neighbor_list>(
                new nl::ContiguousCLLNeighborList(getParticleData(), nl.cll_radius, *pimpl->context, *config)
        );
        log::debug("Using CPU/ContiguousCLL neighbor list with radius {}.", nl.cll_radius);
    } else if (nl.type == "DynamicCLL") {
        (*pimpl).neighborList = std::unique_ptr<neighbor_list>(
                new nl::DynamicCLLNeighborList(getParticleData(), *pimpl->context, *config)
        );
        log::debug("Using CPU/DynamicCLL neighbor list.");
    } else if (nl.type == "Adaptive") {
        (*pimpl).neighborList = std::unique_ptr<neighbor_list>(
                new nl::AdaptiveNeighborList(getParticleData(), *pimpl->context, *config)
        );
        log::debug("Using CPU/Adaptive neighbor list.");
    } else if (nl.type == "CompactCLL") {
        (*pimpl).neighborList = std::unique_ptr<neighbor_list>(
                new nl::CompactCLLNeighborList(getParticleData(), nl.cll_radius, *pimpl->context, *config)
        );
        log::debug("Using CPU/CompactCLL neighbor list with radius {}.", nl.cll_radius);
    } else {
        throw std::invalid_argument(fmt::format(
                R"(In configuration /CPU/NeighborList only "{}", "{}", "{}", or "{}" are valid but not {}.)",
                "CellDecomposition", "ContiguousCLL", "DynamicCLL", "Adaptive", nl.type));
    }
}

scalar &CPUStateModel::energy() {
    return pimpl->currentEnergy;
}

scalar CPUStateModel::energy() const {
    return pimpl->currentEnergy;
}

CPUStateModel::~CPUStateModel() = default;


}
}
}