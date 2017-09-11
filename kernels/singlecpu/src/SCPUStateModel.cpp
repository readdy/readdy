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
 * @file SingleCPUKernelStateModel.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 19.04.16
 */

#include <algorithm>
#include <readdy/kernel/singlecpu/SCPUStateModel.h>

namespace readdy {
namespace kernel {
namespace scpu {
struct SCPUStateModel::Impl {
    using reaction_counts_order1_map = SCPUStateModel::reaction_counts_order1_map;
    using reaction_counts_order2_map = SCPUStateModel::reaction_counts_order2_map;
    scalar currentEnergy = 0;
    model::SCPUParticleData particleData {};
    std::unique_ptr<model::SCPUNeighborList> neighborList;
    SCPUStateModel::topology_action_factory const *topologyActionFactory {nullptr};
    readdy::model::KernelContext const *context {nullptr};
    // only filled when readdy::model::KernelContext::recordReactionsWithPositions is true
    std::vector<readdy::model::reactions::ReactionRecord> reactionRecords{};
    // reaction counts map from particle-type (or type-pair) to vector of count numbers,
    // the position in the vector corresponds to the reaction index,
    // i.e. for each particle-type (or type-pair) there is a new reaction index space
    std::pair<reaction_counts_order1_map, reaction_counts_order2_map> reactionCounts;
};

SCPUStateModel::SCPUStateModel(readdy::model::KernelContext const *context, topology_action_factory const *const taf)
        : pimpl(std::make_unique<SCPUStateModel::Impl>()) {
    pimpl->neighborList = std::make_unique<model::SCPUNeighborList>(context);
    pimpl->context = context;
    pimpl->topologyActionFactory = taf;
}

void readdy::kernel::scpu::SCPUStateModel::addParticle(const readdy::model::Particle &p) {
    pimpl->particleData.addParticles({p});
}

void readdy::kernel::scpu::SCPUStateModel::addParticles(const std::vector<readdy::model::Particle> &p) {
    pimpl->particleData.addParticles(p);
}

const std::vector<readdy::model::Vec3>
readdy::kernel::scpu::SCPUStateModel::getParticlePositions() const {
    const auto &data = pimpl->particleData;
    std::vector<readdy::model::Vec3> target{};
    target.reserve(data.size());
    for (const auto &entry : data) {
        if (!entry.is_deactivated()) target.push_back(entry.position());
    }
    return target;
}

model::SCPUParticleData *readdy::kernel::scpu::SCPUStateModel::getParticleData() const {
    return &pimpl->particleData;
}

void readdy::kernel::scpu::SCPUStateModel::removeParticle(const readdy::model::Particle &p) {
    pimpl->particleData.removeParticle(p);
}

scalar readdy::kernel::scpu::SCPUStateModel::getEnergy() const {
    return pimpl->currentEnergy;
}

void readdy::kernel::scpu::SCPUStateModel::increaseEnergy(scalar increase) {
    pimpl->currentEnergy += increase;
}

const model::SCPUNeighborList *SCPUStateModel::getNeighborList() const {
    return pimpl->neighborList.get();
}

const std::vector<readdy::model::Particle> SCPUStateModel::getParticles() const {
    const auto &data = pimpl->particleData;
    std::vector<readdy::model::Particle> result;
    result.reserve(data.size());
    for (const auto &entry : data) {
        if (!entry.is_deactivated()) {
            result.push_back(data.toParticle(entry));
        }
    }
    return result;
}

void SCPUStateModel::updateNeighborList() {
    pimpl->neighborList->create(pimpl->particleData, pimpl->neighborList->skin());
}

void SCPUStateModel::calculateForces() {
    pimpl->currentEnergy = 0;
    // update forces and energy order 1 potentials
    {
        const readdy::model::Vec3 zero{0, 0, 0};
        for (auto &e : pimpl->particleData) {
            e.force = zero;
            for (const auto &po1 : pimpl->context->potentials().potentials_of(e.type)) {
                po1->calculateForceAndEnergy(e.force, pimpl->currentEnergy, e.position());
            }
        }
    }

    // update forces and energy order 2 potentials
    if(!pimpl->context->potentials().potentials_order2().empty()) {
        const auto &difference = pimpl->context->shortestDifferenceFun();
        readdy::model::Vec3 forceVec{0, 0, 0};
        for (auto it = pimpl->neighborList->begin(); it != pimpl->neighborList->end(); ++it) {
            auto i = it->idx1;
            auto j = it->idx2;
            auto &entry_i = pimpl->particleData.entry_at(i);
            auto &entry_j = pimpl->particleData.entry_at(j);
            const auto &potentials = pimpl->context->potentials().potentials_of(entry_i.type, entry_j.type);
            for (const auto &potential : potentials) {
                potential->calculateForceAndEnergy(forceVec, pimpl->currentEnergy, difference(entry_i.position(), entry_j.position()));
                entry_i.force += forceVec;
                entry_j.force += -1 * forceVec;
            }
        }
    }
    // update forces and energy for topologies
    {
        for ( auto &topology : _topologies) {
            if(!topology->isDeactivated()) {
                // calculate bonded potentials
                for (const auto &bondedPot : topology->getBondedPotentials()) {
                    auto energy = bondedPot->createForceAndEnergyAction(pimpl->topologyActionFactory)->perform(topology.get());
                    pimpl->currentEnergy += energy;
                }
                for (const auto &anglePot : topology->getAnglePotentials()) {
                    auto energy = anglePot->createForceAndEnergyAction(pimpl->topologyActionFactory)->perform(topology.get());
                    pimpl->currentEnergy += energy;
                }
                for (const auto &torsionPot : topology->getTorsionPotentials()) {
                    auto energy = torsionPot->createForceAndEnergyAction(pimpl->topologyActionFactory)->perform(topology.get());
                    pimpl->currentEnergy += energy;
                }
            }
        }
    }
}

void SCPUStateModel::clearNeighborList() {
    pimpl->neighborList->clear();
}

void SCPUStateModel::removeAllParticles() {
    pimpl->particleData.clear();
}

readdy::model::top::GraphTopology *const SCPUStateModel::addTopology(topology_type_type type, const std::vector<readdy::model::TopologyParticle> &particles) {
    std::vector<std::size_t> ids = pimpl->particleData.addTopologyParticles(particles);
    std::vector<particle_type_type> types;
    types.reserve(ids.size());
    for (const auto &p : particles) {
        types.push_back(p.getType());
    }
    auto it = _topologies.emplace_back(
            std::make_unique<topology>(type, std::move(ids), std::move(types),
                                       pimpl->context->topology_registry().potential_configuration())
    );
    const auto idx = std::distance(topologies().begin(), it);
    for(const auto p : (*it)->getParticles()) {
        pimpl->particleData.entry_at(p).topology_index = idx;
    }
    return it->get();
}

std::vector<readdy::model::reactions::ReactionRecord> &SCPUStateModel::reactionRecords() {
    return pimpl->reactionRecords;
}

const std::vector<readdy::model::reactions::ReactionRecord> &SCPUStateModel::reactionRecords() const {
    return pimpl->reactionRecords;
}

std::pair<SCPUStateModel::reaction_counts_order1_map, SCPUStateModel::reaction_counts_order2_map> &SCPUStateModel::reactionCounts() {
    return pimpl->reactionCounts;
}

const std::pair<SCPUStateModel::reaction_counts_order1_map, SCPUStateModel::reaction_counts_order2_map> &SCPUStateModel::reactionCounts() const {
    return pimpl->reactionCounts;
}

readdy::model::Particle SCPUStateModel::getParticleForIndex(const std::size_t index) const {
    return pimpl->particleData.getParticle(index);
}

void SCPUStateModel::expected_n_particles(const std::size_t n) {
    if (pimpl->particleData.size() < n) {
        pimpl->particleData.reserve(n);
    }
}

const SCPUStateModel::topologies_vec &SCPUStateModel::topologies() const {
    return _topologies;
}

SCPUStateModel::topologies_vec &SCPUStateModel::topologies() {
    return _topologies;
}

std::vector<readdy::model::top::GraphTopology*> SCPUStateModel::getTopologies() {
    std::vector<readdy::model::top::GraphTopology*> result;
    result.reserve(_topologies.size() - _topologies.n_deactivated());
    for(const auto& top : _topologies) {
        if(!top->isDeactivated()) {
            result.push_back(top.get());
        }
    }
    return result;
}

particle_type_type SCPUStateModel::getParticleType(const std::size_t index) const {
    return getParticleData()->entry_at(index).type;
}

const readdy::model::top::GraphTopology *SCPUStateModel::getTopologyForParticle(readdy::model::top::Topology::particle_index particle) const {
    const auto& entry = pimpl->particleData.entry_at(particle);
    if(!entry.deactivated) {
        if(entry.topology_index >= 0) {
            return _topologies.at(static_cast<topologies_vec::size_type>(entry.topology_index)).get();
        }
        log::trace("requested particle {} of type {} had no assigned topology", particle, entry.type);
        return nullptr;
    }
    throw std::logic_error(fmt::format("requested particle was deactivated in getTopologyForParticle(p={})", particle));
}

readdy::model::top::GraphTopology *SCPUStateModel::getTopologyForParticle(readdy::model::top::Topology::particle_index particle) {
    const auto& entry = pimpl->particleData.entry_at(particle);
    if(!entry.deactivated) {
        if(entry.topology_index >= 0) {
            return _topologies.at(static_cast<topologies_vec::size_type>(entry.topology_index)).get();
        }
        log::trace("requested particle {} of type {} had no assigned topology", particle, entry.type);
        return nullptr;
    }
    throw std::logic_error(fmt::format("requested particle was deactivated in getTopologyForParticle(p={})", particle));
}

void SCPUStateModel::insert_topology(SCPUStateModel::topology &&top) {
    auto it = _topologies.push_back(std::make_unique<topology>(std::move(top)));
    auto idx = std::distance(_topologies.begin(), it);
    const auto& particles = it->get()->getParticles();
    auto& data = pimpl->particleData;
    std::for_each(particles.begin(), particles.end(), [idx, &data](const topology::particle_index p) {
        data.entry_at(p).topology_index = idx;
    });
}

void SCPUStateModel::initializeNeighborList(scalar skin) {
    pimpl->neighborList->create(pimpl->particleData, skin);
}

SCPUStateModel &SCPUStateModel::operator=(SCPUStateModel &&rhs) noexcept = default;

SCPUStateModel::SCPUStateModel(SCPUStateModel &&rhs) noexcept = default;

SCPUStateModel::~SCPUStateModel() = default;

}
}
}




