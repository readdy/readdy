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
 * @file SingleCPUKernelStateModel.h
 * @brief << brief description >>
 * @author clonker
 * @date 19.04.16
 */

#pragma once

#include <memory>
#include <readdy/model/StateModel.h>
#include <readdy/kernel/singlecpu/model/SCPUParticleData.h>
#include <readdy/model/Context.h>
#include <readdy/kernel/singlecpu/model/SCPUNeighborList.h>
#include <readdy/model/reactions/ReactionRecord.h>
#include "model/ObservableData.h"
#include <readdy/common/index_persistent_vector.h>

namespace readdy {
namespace kernel {
namespace scpu {

class SCPUStateModel : public readdy::model::StateModel {
    using topology_action_factory = readdy::model::top::TopologyActionFactory;
public:
    using reaction_counts_map = readdy::model::reactions::reaction_counts_map;
    using topology = readdy::model::top::GraphTopology;
    using topology_ref = std::unique_ptr<topology>;
    using topologies_vec = readdy::util::index_persistent_vector<topology_ref>;

    void initializeNeighborList(scalar skin) override {
        neighborList->setUp(skin, 1, {});
    }

    void updateNeighborList() override {
        neighborList->update({});
    }

    void clearNeighborList() override {
        neighborList->clear();
    }

    void addParticle(const readdy::model::Particle &p) override {
        particleData.addParticles({p});
    }

    void addParticles(const std::vector<readdy::model::Particle> &p) override {
        particleData.addParticles(p);
    }

    void removeParticle(const readdy::model::Particle &p) override {
        particleData.removeParticle(p);
    }

    void removeAllParticles() override {
        particleData.clear();
    }

    readdy::model::top::GraphTopology *const addTopology(topology_type_type type, const std::vector<readdy::model::TopologyParticle> &particles) override;

    const std::vector<Vec3> getParticlePositions() const override;

    readdy::model::Particle getParticleForIndex(std::size_t index) const override {
        return particleData.getParticle(index);
    }

    particle_type_type getParticleType(std::size_t index) const override {
        return getParticleData()->entry_at(index).type;
    }

    scalar energy() const override {
        return _observableData.energy;
    }

    scalar &energy() override {
        return _observableData.energy;
    }

    SCPUStateModel(const readdy::model::Context &context, const topology_action_factory *);

    ~SCPUStateModel() override = default;

    // move
    SCPUStateModel(SCPUStateModel &&rhs) noexcept = default;

    SCPUStateModel &operator=(SCPUStateModel &&rhs) noexcept = default;

    SCPUStateModel(const SCPUStateModel&) = delete;

    SCPUStateModel& operator=(const SCPUStateModel&) = delete;

    readdy::kernel::scpu::model::SCPUParticleData *getParticleData() {
        return &particleData;
    }

    const readdy::kernel::scpu::model::SCPUParticleData *getParticleData() const {
        return &particleData;
    }

    virtual const model::CellLinkedList *getNeighborList() const {
        return neighborList.get();
    }

    model::CellLinkedList *getNeighborList() {
        return neighborList.get();
    }

    const std::vector<readdy::model::Particle> getParticles() const override;

    std::vector<readdy::model::reactions::ReactionRecord>& reactionRecords() {
        return _observableData.reactionRecords;
    }

    const std::vector<readdy::model::reactions::ReactionRecord>& reactionRecords() const {
        return _observableData.reactionRecords;
    }

    const reaction_counts_map & reactionCounts() const {
        return _observableData.reactionCounts;
    }

    reaction_counts_map &reactionCounts() {
        return _observableData.reactionCounts;
    }

    Matrix33 &virial() {
        return _observableData.virial;
    }

    const Matrix33 &virial() const {
        return _observableData.virial;
    }

    void resetReactionCounts();

    const topologies_vec &topologies() const {
        return _topologies;
    }

    topologies_vec &topologies() {
        return _topologies;
    }

    model::ObservableData &observableData() {
        return _observableData;
    }

    const model::ObservableData &observableData() const {
        return _observableData;
    }

    void insert_topology(topology&& top);

    void toDenseParticleIndices(std::vector<std::size_t>::iterator begin,
                                std::vector<std::size_t>::iterator end) const override;

    std::vector<readdy::model::top::GraphTopology*> getTopologies() override;

    const readdy::model::top::GraphTopology *getTopologyForParticle(readdy::model::top::Topology::particle_index particle) const override;

    readdy::model::top::GraphTopology *getTopologyForParticle(readdy::model::top::Topology::particle_index particle) override;

private:
    model::SCPUParticleData particleData {};
    std::unique_ptr<model::CellLinkedList> neighborList;
    SCPUStateModel::topology_action_factory const *topologyActionFactory {nullptr};
    // only filled when readdy::model::Context::recordReactionsWithPositions is true

    std::reference_wrapper<const readdy::model::Context> _context;

    topologies_vec _topologies;

    model::ObservableData _observableData;
};

}
}
}
