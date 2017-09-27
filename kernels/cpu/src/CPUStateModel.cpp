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
using dist_fun = readdy::model::Context::shortest_dist_fun;

const std::vector<Vec3> CPUStateModel::getParticlePositions() const {
    const auto data = _neighborList->data();
    std::vector<Vec3> target{};
    target.reserve(data->size());
    for (const auto &entry : *data) {
        if (!entry.deactivated) target.push_back(entry.pos);
    }
    return target;
}

const std::vector<readdy::model::Particle> CPUStateModel::getParticles() const {
    const auto data = _neighborList->data();
    std::vector<readdy::model::Particle> result;
    result.reserve(data->size());
    for (const auto &entry : *data) {
        if (!entry.deactivated) {
            result.push_back(data->toParticle(entry));
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

CPUStateModel::CPUStateModel(const readdy::model::Context &context,
                             readdy::util::thread::Config const *const config,
                             readdy::model::top::TopologyActionFactory const *const taf)
        : _config(*config), _context(context), _topologyActionFactory(*taf) {
    _neighborList = std::unique_ptr<neighbor_list>(
            // new nl::CellDecompositionNeighborList(*getParticleData(), *pimpl->context, *config)
            new nl::CompactCLLNeighborList(1, _context.get(), *config)
            //new nl::DynamicCLLNeighborList(*pimpl->context, *config)
    );
    _reorderConnection = std::make_unique<readdy::signals::scoped_connection>(
            getParticleData()->registerReorderEventListener([this](const std::vector<std::size_t> &indices) -> void {
                for (auto &top : _topologies) {
                    if(!top->isDeactivated()) top->permuteIndices(indices);
                }
            }));
}

CPUStateModel::data_type const *const CPUStateModel::getParticleData() const {
    return _neighborList->data();
}

CPUStateModel::data_type *const CPUStateModel::getParticleData() {
    return _neighborList->data();
}

neighbor_list const *const CPUStateModel::getNeighborList() const {
    return _neighborList.get();
}

void CPUStateModel::clearNeighborList() {
    clearNeighborList({});
}

void CPUStateModel::removeAllParticles() {
    getParticleData()->clear();
}

neighbor_list *const CPUStateModel::getNeighborList() {
    return _neighborList.get();
}

readdy::model::top::GraphTopology *const
CPUStateModel::addTopology(topology_type_type type, const std::vector<readdy::model::TopologyParticle> &particles) {
    std::vector<std::size_t> ids = getParticleData()->addTopologyParticles(particles);
    std::vector<particle_type_type> types;
    types.reserve(ids.size());
    for (const auto &p : particles) {
        types.push_back(p.getType());
    }
    auto it = _topologies.push_back(std::make_unique<topology>(type, std::move(ids), std::move(types), _context.get()));
    const auto idx = std::distance(topologies().begin(), it);
    for(const auto p : (*it)->getParticles()) {
        getParticleData()->entry_at(p).topology_index = idx;
    }
    return it->get();
}

std::vector<readdy::model::reactions::ReactionRecord> &CPUStateModel::reactionRecords() {
    return _reactionRecords;
}

const std::vector<readdy::model::reactions::ReactionRecord> &CPUStateModel::reactionRecords() const {
    return _reactionRecords;
}

readdy::model::Particle CPUStateModel::getParticleForIndex(const std::size_t index) const {
    return getParticleData()->getParticle(index);
}

const std::pair<CPUStateModel::reaction_counts_order1_map, CPUStateModel::reaction_counts_order2_map> &CPUStateModel::reactionCounts() const {
    return _reactionCounts;
}

std::pair<CPUStateModel::reaction_counts_order1_map, CPUStateModel::reaction_counts_order2_map> &CPUStateModel::reactionCounts() {
    return _reactionCounts;
}

const CPUStateModel::topologies_vec &CPUStateModel::topologies() const {
    return _topologies;
}

CPUStateModel::topologies_vec &CPUStateModel::topologies() {
    return _topologies;
}

std::vector<readdy::model::top::GraphTopology*> CPUStateModel::getTopologies() {
    std::vector<readdy::model::top::GraphTopology*> result;
    result.reserve(_topologies.size() - _topologies.n_deactivated());
    for(const auto& top : _topologies) {
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
            return _topologies.at(static_cast<topologies_vec::size_type>(entry.topology_index)).get();
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
            return _topologies.at(static_cast<topologies_vec::size_type>(entry.topology_index)).get();
        }
        log::trace("requested particle {} of type {} had no assigned topology", particle, entry.type);
        return nullptr;
    }
    throw std::logic_error(fmt::format("requested particle was deactivated in getTopologyForParticle(p={})", particle));
}

void CPUStateModel::insert_topology(CPUStateModel::topology &&top) {
    auto it = _topologies.push_back(std::make_unique<topology>(std::move(top)));
    auto idx = std::distance(_topologies.begin(), it);
    const auto& particles = it->get()->getParticles();
    auto& data = *getParticleData();
    std::for_each(particles.begin(), particles.end(), [idx, &data](const topology::particle_index p) {
        data.entry_at(p).topology_index = idx;
    });
}

void CPUStateModel::updateNeighborList(const util::PerformanceNode &node) {
    _neighborList->update(node.subnode("update"));
}

void CPUStateModel::clearNeighborList(const util::PerformanceNode &node) {
    _neighborList->clear(node.subnode("clear"));
}

void CPUStateModel::initializeNeighborList(scalar skin) {
    initializeNeighborList(skin, {});
}

void CPUStateModel::initializeNeighborList(scalar skin, const util::PerformanceNode &node) {
    _neighborList->skin() = skin;
    _neighborList->set_up(node.subnode("set_up"));
}

void CPUStateModel::configure(const readdy::conf::cpu::Configuration &configuration) {
    const auto& nl = configuration.neighborList;

    if(nl.type == "CellDecomposition") {
        _neighborList = std::unique_ptr<neighbor_list>(
                new nl::CellDecompositionNeighborList(getParticleData(), _context.get(), _config.get())
        );
        log::debug("Using CPU/CellDecomposition neighborList");
    } else if(nl.type == "ContiguousCLL") {
        _neighborList = std::unique_ptr<neighbor_list>(
                new nl::ContiguousCLLNeighborList(getParticleData(), nl.cll_radius, _context.get(), _config.get())
        );
        log::debug("Using CPU/ContiguousCLL neighbor list with radius {}.", nl.cll_radius);
    } else if (nl.type == "DynamicCLL") {
        _neighborList = std::unique_ptr<neighbor_list>(
                new nl::DynamicCLLNeighborList(getParticleData(), _context.get(), _config.get())
        );
        log::debug("Using CPU/DynamicCLL neighbor list.");
    } else if (nl.type == "Adaptive") {
        _neighborList = std::unique_ptr<neighbor_list>(
                new nl::AdaptiveNeighborList(getParticleData(), _context.get(), _config.get())
        );
        log::debug("Using CPU/Adaptive neighbor list.");
    } else if (nl.type == "CompactCLL") {
        _neighborList = std::unique_ptr<neighbor_list>(
                new nl::CompactCLLNeighborList(getParticleData(), nl.cll_radius, _context.get(), _config.get())
        );
        log::debug("Using CPU/CompactCLL neighbor list with radius {}.", nl.cll_radius);
    } else {
        throw std::invalid_argument(fmt::format(
                R"(In configuration /CPU/NeighborList only "{}", "{}", "{}", or "{}" are valid but not {}.)",
                "CellDecomposition", "ContiguousCLL", "DynamicCLL", "Adaptive", nl.type));
    }
}

scalar &CPUStateModel::energy() {
    return _currentEnergy;
}

scalar CPUStateModel::energy() const {
    return _currentEnergy;
}

CPUStateModel::~CPUStateModel() = default;


}
}
}