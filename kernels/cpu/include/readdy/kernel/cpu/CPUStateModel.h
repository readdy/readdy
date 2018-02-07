/********************************************************************
 * Copyright © 2017 Computational Molecular Biology Group,          *
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
 * @file CPUStateModel.h
 * @brief << brief description >>
 * @author clonker
 * @date 12/11/17
 */



#pragma once

#include <readdy/model/StateModel.h>


#include <readdy/model/StateModel.h>
#include <readdy/model/Context.h>
#include <readdy/common/thread/Config.h>
#include <readdy/model/reactions/ReactionRecord.h>
#include <readdy/model/observables/ReactionCounts.h>
#include <readdy/common/index_persistent_vector.h>
#include <readdy/common/Timer.h>
#include <readdy/api/KernelConfiguration.h>
#include <readdy/kernel/cpu/data/DefaultDataContainer.h>
#include <readdy/kernel/cpu/nl/CellLinkedList.h>
#include <readdy/kernel/cpu/nl/ContiguousCellLinkedList.h>
#include <readdy/kernel/cpu/data/ObservableData.h>

namespace readdy {
namespace kernel {
namespace cpu {
class CPUStateModel : public readdy::model::StateModel {

public:

    using data_type = readdy::kernel::cpu::data::DefaultDataContainer;
    using particle_type = readdy::model::Particle;
    using reaction_counts_map = readdy::model::reactions::reaction_counts_map;

    using topology = readdy::model::top::GraphTopology;
    using topology_ref = std::unique_ptr<topology>;
    using topologies_vec = readdy::util::index_persistent_vector<topology_ref>;
    using neighbor_list = nl::CompactCellLinkedList;

    CPUStateModel(data_type &data, const readdy::model::Context &context, thread_pool &pool,
                  readdy::model::top::TopologyActionFactory const* taf);

    ~CPUStateModel() override = default;

    CPUStateModel(const CPUStateModel&) = delete;
    CPUStateModel& operator=(const CPUStateModel&) = delete;
    CPUStateModel(CPUStateModel&&) = delete;
    CPUStateModel& operator=(CPUStateModel&&) = delete;

    void configure(const readdy::conf::cpu::Configuration &configuration) {
        const auto& nl = configuration.neighborList;
        _neighborListCellRadius = nl.cll_radius;
    }

    const std::vector<Vec3> getParticlePositions() const override;

    const std::vector<particle_type> getParticles() const override;

    void initializeNeighborList(scalar skin, const util::PerformanceNode &node) {
        _neighborList->setUp(skin, _neighborListCellRadius, node.subnode("set_up"));
        _neighborList->update(node.subnode("update"));
    };

    void initializeNeighborList(scalar skin) override {
        initializeNeighborList(skin, {});
    };

    void updateNeighborList(const util::PerformanceNode &node) {
        _neighborList->update(node.subnode("update"));
    };

    void updateNeighborList() override {
        updateNeighborList({});
    };

    void addParticle(const particle_type &p) override {
        getParticleData()->addParticle(p);
    };

    void addParticles(const std::vector<particle_type> &p) override {
        getParticleData()->addParticles(p);
    };

    void removeParticle(const particle_type &p) override {
        getParticleData()->removeParticle(p);
    };

    void removeAllParticles() override {
        getParticleData()->clear();
    };

    data::ObservableData &observableData() {
        return _observableData;
    }
    
    const data::ObservableData &observableData() const {
        return _observableData;
    }
    
    Matrix33 &virial() {
        return _observableData.virial;
    }
    
    const Matrix33 &virial() const {
        return _observableData.virial;
    }
    
    scalar energy() const override {
        return _observableData.energy;
    };

    scalar &energy() override {
        return _observableData.energy;
    };

    data_type const *const getParticleData() const {
        return &_data.get();
    };

    data_type *const getParticleData() {
        return &_data.get();
    };

    neighbor_list const *const getNeighborList() const {
        return _neighborList.get();

    };

    neighbor_list *const getNeighborList() {
        return _neighborList.get();
    };

    void clearNeighborList() override {
        clearNeighborList({});
    };

    void clearNeighborList(const util::PerformanceNode &node) {
        _neighborList->clear();
    };

    readdy::model::top::GraphTopology *const
    addTopology(topology_type_type type, const std::vector<readdy::model::TopologyParticle> &particles) override;

    std::vector<readdy::model::reactions::ReactionRecord> &reactionRecords() {
        return _observableData.reactionRecords;
    };

    const std::vector<readdy::model::reactions::ReactionRecord> &reactionRecords() const {
        return _observableData.reactionRecords;
    };

    const reaction_counts_map & reactionCounts() const {
        return _observableData.reactionCounts;
    };

    reaction_counts_map &reactionCounts() {
        return _observableData.reactionCounts;
    };

    void resetReactionCounts();

    particle_type getParticleForIndex(std::size_t index) const override {
        return _data.get().getParticle(index);
    };

    particle_type_type getParticleType(std::size_t index) const override {
        return _data.get().entry_at(index).type;
    };

    const topologies_vec &topologies() const {
        return _topologies;
    };

    topologies_vec &topologies() {
        return _topologies;
    };

    void insert_topology(topology&& top);

    std::vector<readdy::model::top::GraphTopology *> getTopologies() override;

    const readdy::model::top::GraphTopology *getTopologyForParticle(readdy::model::top::Topology::particle_index particle) const override;

    readdy::model::top::GraphTopology *getTopologyForParticle(readdy::model::top::Topology::particle_index particle) override;

    void toDenseParticleIndices(std::vector<std::size_t>::iterator begin,
                                std::vector<std::size_t>::iterator end) const override;

private:
    data::ObservableData _observableData;
    std::reference_wrapper<thread_pool> _pool;
    std::reference_wrapper<const readdy::model::Context> _context;
    std::reference_wrapper<data_type> _data;
    std::unique_ptr<neighbor_list> _neighborList;
    neighbor_list::cell_radius_type _neighborListCellRadius {1};
    std::unique_ptr<readdy::signals::scoped_connection> _reorderConnection;
    std::reference_wrapper<const readdy::model::top::TopologyActionFactory> _topologyActionFactory;
    topologies_vec _topologies{};
};
}
}
}
