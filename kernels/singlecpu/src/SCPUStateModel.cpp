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
    using reaction_counts_map = SCPUStateModel::reaction_counts_map;
    scalar currentEnergy = 0;
    model::SCPUParticleData particleData {};
    std::unique_ptr<model::SCPUNeighborList> neighborList;
    SCPUStateModel::topology_action_factory const *topologyActionFactory {nullptr};
    // only filled when readdy::model::Context::recordReactionsWithPositions is true
    std::vector<readdy::model::reactions::ReactionRecord> reactionRecords{};
    reaction_counts_map reactionCounts {};
};

SCPUStateModel::SCPUStateModel(const readdy::model::Context &context, topology_action_factory const *const taf)
        : pimpl(std::make_unique<SCPUStateModel::Impl>()), _context(context) {
    pimpl->neighborList = std::make_unique<model::SCPUNeighborList>(&_context.get());
    pimpl->topologyActionFactory = taf;
}

void readdy::kernel::scpu::SCPUStateModel::addParticle(const readdy::model::Particle &p) {
    pimpl->particleData.addParticles({p});
}

void readdy::kernel::scpu::SCPUStateModel::addParticles(const std::vector<readdy::model::Particle> &p) {
    pimpl->particleData.addParticles(p);
}

const std::vector<Vec3>
readdy::kernel::scpu::SCPUStateModel::getParticlePositions() const {
    const auto &data = pimpl->particleData;
    std::vector<Vec3> target{};
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
            std::make_unique<topology>(type, std::move(ids), std::move(types), _context.get(), this)
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

SCPUStateModel::reaction_counts_map &SCPUStateModel::reactionCounts() {
    return pimpl->reactionCounts;
}

const SCPUStateModel::reaction_counts_map &SCPUStateModel::reactionCounts() const {
    return pimpl->reactionCounts;
}

readdy::model::Particle SCPUStateModel::getParticleForIndex(const std::size_t index) const {
    return pimpl->particleData.getParticle(index);
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

scalar SCPUStateModel::energy() const {
    return pimpl->currentEnergy;
}

scalar &SCPUStateModel::energy() {
    return pimpl->currentEnergy;
}

SCPUStateModel &SCPUStateModel::operator=(SCPUStateModel &&rhs) noexcept = default;

SCPUStateModel::SCPUStateModel(SCPUStateModel &&rhs) noexcept = default;

SCPUStateModel::~SCPUStateModel() = default;

}
}
}




