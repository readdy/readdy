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
 * << detailed description >>
 *
 * @file SCPUStateModel.h
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

    void initializeNeighborList(scalar interactionDistance) override {
        neighborList->setUp(interactionDistance, 1);
    }

    void updateNeighborList() override {
        neighborList->update();
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

    readdy::model::top::GraphTopology *const addTopology(TopologyTypeId type, const std::vector<readdy::model::TopologyParticle> &particles) override;

    const std::vector<Vec3> getParticlePositions() const override;

    readdy::model::Particle getParticleForIndex(std::size_t index) const override {
        return particleData.getParticle(index);
    }

    ParticleTypeId getParticleType(std::size_t index) const override {
        return getParticleData()->entry_at(index).type;
    }

    scalar energy() const override {
        return _observableData.energy;
    }

    scalar &energy() override {
        return _observableData.energy;
    }

    scalar time() const override {
        return _observableData.time;
    }

    scalar &time() override {
        return _observableData.time;
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

    void clear() override;

private:
    model::SCPUParticleData particleData {};
    std::unique_ptr<model::CellLinkedList> neighborList;
    SCPUStateModel::topology_action_factory const *topologyActionFactory {nullptr};

    std::reference_wrapper<const readdy::model::Context> _context;

    topologies_vec _topologies;

    model::ObservableData _observableData;
};

}
}
}
