/********************************************************************
 * Copyright © 2019 Computational Molecular Biology Group,          *
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
 * « detailed description »
 *
 * @file MPIStateModel.h
 * @brief « brief description »
 * @author chrisfroe
 * @date 28.05.19
 */

#pragma once

#include <readdy/model/StateModel.h>
#include <readdy/kernel/mpi/model/MPINeighborList.h>
#include <readdy/kernel/singlecpu/model/ObservableData.h>
#include <readdy/model/reactions/ReactionRecord.h>
#include <readdy/common/signals.h>
#include <readdy/kernel/mpi/model/MPIParticleData.h>

namespace readdy::kernel::mpi {

class MPIStateModel : public readdy::model::StateModel {

public:

    using Data = MPIDataContainer;
    using Particle = readdy::model::Particle;
    using ReactionCountsMap = readdy::model::reactions::reaction_counts_map;
    using NeighborList = model::CellLinkedList;

    MPIStateModel(Data &data, const readdy::model::Context &context);

    ~MPIStateModel() override = default;

    MPIStateModel(const MPIStateModel&) = delete;
    MPIStateModel& operator=(const MPIStateModel&) = delete;
    MPIStateModel(MPIStateModel&&) = delete;
    MPIStateModel& operator=(MPIStateModel&&) = delete;

    const std::vector<Vec3> getParticlePositions() const override;

    const std::vector<Particle> getParticles() const override;

    void initializeNeighborList(scalar interactionDistance) override {
        _neighborList->setUp(interactionDistance, _neighborListCellRadius);
        _neighborList->update();
    };

    void updateNeighborList() override {
        _neighborList->update();
    };

    void addParticle(const Particle &p) override {
        getParticleData()->addParticle(p);
    };

    void addParticles(const std::vector<Particle> &p) override {
        getParticleData()->addParticles(p);
    };

    void removeParticle(const Particle &p) override {
        getParticleData()->removeParticle(p);
    };

    void removeAllParticles() override {
        getParticleData()->clear();
    };

    readdy::kernel::scpu::model::ObservableData &observableData() {
        return _observableData;
    }

    const readdy::kernel::scpu::model::ObservableData &observableData() const {
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

    scalar &time() override {
        return _observableData.time;
    };

    Data const *const getParticleData() const {
        return &_data.get();
    };

    Data *const getParticleData() {
        return &_data.get();
    };

    NeighborList const *const getNeighborList() const {
        return _neighborList.get();

    };

    NeighborList *const getNeighborList() {
        return _neighborList.get();
    };

    void clearNeighborList() override {
        _neighborList->clear();
    };

    std::vector<readdy::model::reactions::ReactionRecord> &reactionRecords() {
        return _observableData.reactionRecords;
    };

    const std::vector<readdy::model::reactions::ReactionRecord> &reactionRecords() const {
        return _observableData.reactionRecords;
    };

    const ReactionCountsMap & reactionCounts() const {
        return _observableData.reactionCounts;
    };

    ReactionCountsMap &reactionCounts() {
        return _observableData.reactionCounts;
    };

    void resetReactionCounts();

    Particle getParticleForIndex(std::size_t index) const override {
        return _data.get().getParticle(index);
    };

    ParticleTypeId getParticleType(std::size_t index) const override {
        return _data.get().entry_at(index).type;
    };

    void toDenseParticleIndices(std::vector<std::size_t>::iterator begin,
                                std::vector<std::size_t>::iterator end) const override;

    void clear() override;

    readdy::model::top::GraphTopology *const
    addTopology(TopologyTypeId type, const std::vector<readdy::model::TopologyParticle> &particles) override {
        throw std::logic_error("no topologies on MPI kernel");
    }

    std::vector<readdy::model::top::GraphTopology *> getTopologies() override {
        throw std::logic_error("no topologies on MPI kernel");
    }

    const readdy::model::top::GraphTopology *getTopologyForParticle(readdy::model::top::Topology::particle_index particle) const override {
        throw std::logic_error("no topologies on MPI kernel");
    }

    readdy::model::top::GraphTopology *getTopologyForParticle(readdy::model::top::Topology::particle_index particle) override {
        throw std::logic_error("no topologies on MPI kernel");
    }

private:
    readdy::kernel::scpu::model::ObservableData _observableData;
    std::reference_wrapper<const readdy::model::Context> _context;
    std::reference_wrapper<Data> _data;
    std::unique_ptr<NeighborList> _neighborList;
    NeighborList::CellRadius _neighborListCellRadius {1};
    std::unique_ptr<readdy::signals::scoped_connection> _reorderConnection;
};

}
