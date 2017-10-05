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
#include <readdy/common/index_persistent_vector.h>
#include <readdy/model/reactions/Utils.h>

namespace readdy {
namespace kernel {
namespace scpu {

class SCPUStateModel : public readdy::model::StateModel {
    using topology_action_factory = readdy::model::top::TopologyActionFactory;
public:
    using reaction_counts_map = readdy::model::reactions::utils::reaction_counts_map;
    using topology = readdy::model::top::GraphTopology;
    using topology_ref = std::unique_ptr<topology>;
    using topologies_vec = readdy::util::index_persistent_vector<topology_ref>;

    void initializeNeighborList(scalar skin) override;

    void updateNeighborList() override;

    void clearNeighborList() override;

    void addParticle(const readdy::model::Particle &p) override;

    void addParticles(const std::vector<readdy::model::Particle> &p) override;

    void removeParticle(const readdy::model::Particle &p) override;

    void removeAllParticles() override;

    readdy::model::top::GraphTopology *const addTopology(topology_type_type type, const std::vector<readdy::model::TopologyParticle> &particles) override;

    const std::vector<Vec3> getParticlePositions() const override;

    readdy::model::Particle getParticleForIndex(std::size_t index) const override;

    particle_type_type getParticleType(std::size_t index) const override;

    scalar energy() const override;

    scalar &energy() override;

    virtual void increaseEnergy(scalar increase);

    SCPUStateModel(const readdy::model::Context &context, const topology_action_factory *);

    ~SCPUStateModel() override;

    // move
    SCPUStateModel(SCPUStateModel &&rhs) noexcept;

    SCPUStateModel &operator=(SCPUStateModel &&rhs) noexcept;

    SCPUStateModel(const SCPUStateModel&) = delete;

    SCPUStateModel& operator=(const SCPUStateModel&) = delete;

    virtual readdy::kernel::scpu::model::SCPUParticleData *getParticleData() const;

    virtual const model::SCPUNeighborList *getNeighborList() const;

    const std::vector<readdy::model::Particle> getParticles() const override;

    std::vector<readdy::model::reactions::ReactionRecord>& reactionRecords();

    const std::vector<readdy::model::reactions::ReactionRecord>& reactionRecords() const;

    const reaction_counts_map & reactionCounts() const;

    reaction_counts_map &reactionCounts();

    const topologies_vec &topologies() const;

    topologies_vec &topologies();

    void insert_topology(topology&& top);

    std::vector<readdy::model::top::GraphTopology*> getTopologies() override;

    const readdy::model::top::GraphTopology *getTopologyForParticle(readdy::model::top::Topology::particle_index particle) const override;

    readdy::model::top::GraphTopology *getTopologyForParticle(readdy::model::top::Topology::particle_index particle) override;

private:
    struct Impl;
    std::unique_ptr<Impl> pimpl;

    std::reference_wrapper<const readdy::model::Context> _context;

    topologies_vec _topologies;
};

}
}
}
