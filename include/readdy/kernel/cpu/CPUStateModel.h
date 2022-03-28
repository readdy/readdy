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
#include <readdy/api/KernelConfiguration.h>
#include <readdy/kernel/cpu/data/DefaultDataContainer.h>
#include <readdy/kernel/cpu/nl/CellLinkedList.h>
#include <readdy/kernel/cpu/nl/ContiguousCellLinkedList.h>
#include <readdy/kernel/cpu/data/ObservableData.h>

namespace readdy::kernel::cpu {
class CPUStateModel : public readdy::model::StateModel {

public:

    using data_type = readdy::kernel::cpu::data::DefaultDataContainer;
    using particle_type = readdy::model::Particle;
    
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

    std::vector<Vec3> getParticlePositions() const override;

    std::vector<particle_type> getParticles() const override;

    void initializeNeighborList(scalar interactionDistance) override {
        _neighborList->setUp(interactionDistance, _neighborListCellRadius);
        _neighborList->update();
    };

    void updateNeighborList() override {
        _neighborList->update();
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

    scalar time() const override {
        return _observableData.time;
    };

    void setTime(scalar t) override {
        _observableData.time = t;
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
        _neighborList->clear();
    };

    readdy::model::top::GraphTopology *const
    addTopology(TopologyTypeId type, const std::vector<readdy::model::Particle> &particles) override;

    std::vector<readdy::model::reactions::ReactionRecord> &reactionRecords() {
        return _observableData.reactionRecords;
    };

    const std::vector<readdy::model::reactions::ReactionRecord> &reactionRecords() const {
        return _observableData.reactionRecords;
    };

    const readdy::model::reactions::ReactionCounts & reactionCounts() const {
        return _observableData.reactionCounts;
    };

    readdy::model::reactions::ReactionCounts &reactionCounts() {
        return _observableData.reactionCounts;
    };

    const readdy::model::reactions::SpatialTopologyReactionCounts &spatialReactionCounts() const {
        return _observableData.spatialReactionCounts;
    }

    readdy::model::reactions::SpatialTopologyReactionCounts &spatialReactionCounts() {
        return _observableData.spatialReactionCounts;
    }

    const readdy::model::reactions::StructuralTopologyReactionCounts &structuralReactionCounts() const {
        return _observableData.structuralReactionCounts;
    }

    readdy::model::reactions::StructuralTopologyReactionCounts &structuralReactionCounts() {
        return _observableData.structuralReactionCounts;
    }

    void resetReactionCounts();
    
    void resetTopologyReactionCounts() {
        const auto &topologies = _context.get().topologyRegistry();
        for (const auto &entry : topologies.spatialReactionRegistry()) {
            for (const auto &sr : entry.second) {
                _observableData.spatialReactionCounts[sr.id()] = 0;
            }
        }
        for(const auto &type : topologies.types()) {
            for(const auto &structuralReaction : type.structuralReactions) {
                _observableData.structuralReactionCounts[structuralReaction.id()] = 0;
            }
        }
    }

    particle_type getParticleForIndex(std::size_t index) const override {
        return _data.get().getParticle(index);
    };

    ParticleTypeId getParticleType(std::size_t index) const override {
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

    void toDenseParticleIndices(std::vector<std::size_t>::iterator begin,
                                std::vector<std::size_t>::iterator end) const override;

    void clear() override;

private:
    data::ObservableData _observableData;
    std::reference_wrapper<thread_pool> _pool;
    std::reference_wrapper<const readdy::model::Context> _context;
    std::reference_wrapper<data_type> _data;
    std::unique_ptr<neighbor_list> _neighborList;
    neighbor_list::cell_radius_type _neighborListCellRadius {1};
    std::reference_wrapper<const readdy::model::top::TopologyActionFactory> _topologyActionFactory;
    topologies_vec _topologies{};
};
}
