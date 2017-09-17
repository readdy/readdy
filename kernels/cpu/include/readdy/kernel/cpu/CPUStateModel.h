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
 * @file CPUStateModel.h
 * @brief << brief description >>
 * @author clonker
 * @date 13.07.16
 */

#pragma once

#include <readdy/model/KernelStateModel.h>
#include <readdy/model/KernelContext.h>
#include <readdy/common/thread/Config.h>
#include <readdy/model/reactions/ReactionRecord.h>
#include <readdy/model/observables/ReactionCounts.h>
#include <readdy/kernel/cpu/util/config.h>
#include <readdy/common/index_persistent_vector.h>
#include <readdy/common/Timer.h>
#include <readdy/api/KernelConfiguration.h>
#include <readdy/kernel/cpu/data/DefaultDataContainer.h>

namespace readdy {
namespace kernel {
namespace cpu {
class CPUStateModel : public readdy::model::KernelStateModel {

public:

    using data_type = readdy::kernel::cpu::data::EntryDataContainer;
    using particle_type = readdy::model::Particle;
    using reaction_counts_order1_map = readdy::model::observables::ReactionCounts::reaction_counts_order1_map;
    using reaction_counts_order2_map = readdy::model::observables::ReactionCounts::reaction_counts_order2_map;

    using topology = readdy::model::top::GraphTopology;
    using topology_ref = std::unique_ptr<topology>;
    using topologies_vec = readdy::util::index_persistent_vector<topology_ref>;

    CPUStateModel(readdy::model::KernelContext* context, readdy::util::thread::Config const* config,
                  readdy::model::top::TopologyActionFactory const* taf);

    ~CPUStateModel() override;

    CPUStateModel(const CPUStateModel&) = delete;
    CPUStateModel& operator=(const CPUStateModel&) = delete;
    CPUStateModel(CPUStateModel&&) = delete;
    CPUStateModel& operator=(CPUStateModel&&) = delete;

    void configure(const readdy::conf::cpu::Configuration &configuration);

    const std::vector<Vec3> getParticlePositions() const override;

    const std::vector<particle_type> getParticles() const override;

    void initializeNeighborList(scalar skin, const util::PerformanceNode &node);

    void initializeNeighborList(scalar skin) override;

    void updateNeighborList(const util::PerformanceNode &node);

    void updateNeighborList() override;

    void calculateForces() override;

    void addParticle(const particle_type &p) override;

    void addParticles(const std::vector<particle_type> &p) override;

    void removeParticle(const particle_type &p) override;

    void removeAllParticles() override;

    readdy::scalar getEnergy() const override;

    data_type const *const getParticleData() const;

    data_type *const getParticleData();

    neighbor_list const *const getNeighborList() const;

    neighbor_list *const getNeighborList();

    void clearNeighborList() override;

    void clearNeighborList(const util::PerformanceNode &node);

    readdy::model::top::GraphTopology *const
    addTopology(topology_type_type type, const std::vector<readdy::model::TopologyParticle> &particles) override;

    std::vector<readdy::model::reactions::ReactionRecord> &reactionRecords();

    const std::vector<readdy::model::reactions::ReactionRecord> &reactionRecords() const;

    const std::pair<reaction_counts_order1_map, reaction_counts_order2_map> &reactionCounts() const;

    std::pair<reaction_counts_order1_map, reaction_counts_order2_map> &reactionCounts();

    particle_type getParticleForIndex(std::size_t index) const override;

    particle_type_type getParticleType(std::size_t index) const override;

    const topologies_vec &topologies() const;

    topologies_vec &topologies();

    void insert_topology(topology&& top);

    std::vector<readdy::model::top::GraphTopology *> getTopologies() override;

    const readdy::model::top::GraphTopology *getTopologyForParticle(readdy::model::top::Topology::particle_index particle) const override;

    readdy::model::top::GraphTopology *getTopologyForParticle(readdy::model::top::Topology::particle_index particle) override;

private:
    struct Impl;
    std::unique_ptr<Impl> pimpl;
    readdy::util::thread::Config const *const config;
};
}
}
}
