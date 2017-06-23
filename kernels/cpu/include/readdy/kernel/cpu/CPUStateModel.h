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
#include <readdy/kernel/cpu/model/ParticleIndexPair.h>
#include <readdy/common/thread/Config.h>
#include <readdy/kernel/cpu/model/CPUParticleData.h>
#include <readdy/model/reactions/ReactionRecord.h>
#include <readdy/model/observables/ReactionCounts.h>
#include <readdy/kernel/cpu/util/config.h>
#include <readdy/common/index_persistent_vector.h>

namespace readdy {
namespace kernel {
namespace cpu {
class CPUStateModel : public readdy::model::KernelStateModel {

public:

    using data_t = readdy::kernel::cpu::model::CPUParticleData;
    using particle_t = readdy::model::Particle;
    using reaction_counts_order1_map = readdy::model::observables::ReactionCounts::reaction_counts_order1_map;
    using reaction_counts_order2_map = readdy::model::observables::ReactionCounts::reaction_counts_order2_map;

    using topologies_t = readdy::util::index_persistent_vector<readdy::model::top::GraphTopology>;

    CPUStateModel(readdy::model::KernelContext* context, readdy::util::thread::Config const* config,
                  readdy::model::top::TopologyActionFactory const* taf);

    ~CPUStateModel() override;

    CPUStateModel(const CPUStateModel&) = delete;
    CPUStateModel& operator=(const CPUStateModel&) = delete;
    CPUStateModel(CPUStateModel&&) = delete;
    CPUStateModel& operator=(CPUStateModel&&) = delete;

    const std::vector<readdy::model::Vec3> getParticlePositions() const override;

    const std::vector<particle_t> getParticles() const override;

    void updateNeighborList() override;

    void calculateForces() override;

    void addParticle(const particle_t &p) override;

    void addParticles(const std::vector<particle_t> &p) override;

    void removeParticle(const particle_t &p) override;

    void removeAllParticles() override;

    readdy::scalar getEnergy() const override;

    data_t const *const getParticleData() const;

    data_t *const getParticleData();

    neighbor_list const *const getNeighborList() const;

    neighbor_list *const getNeighborList();

    void expected_n_particles(std::size_t n) override;

    void clearNeighborList() override;

    readdy::model::top::GraphTopology *const
    addTopology(const std::vector<readdy::model::TopologyParticle> &particles) override;

    std::vector<readdy::model::reactions::ReactionRecord> &reactionRecords();

    const std::vector<readdy::model::reactions::ReactionRecord> &reactionRecords() const;

    const std::pair<reaction_counts_order1_map, reaction_counts_order2_map> &reactionCounts() const;

    std::pair<reaction_counts_order1_map, reaction_counts_order2_map> &reactionCounts();

    particle_t getParticleForIndex(std::size_t index) const override;

    particle_type_type getParticleType(std::size_t index) const override;

    const topologies_t &topologies() const;

    topologies_t &topologies();

    virtual std::vector<readdy::model::top::GraphTopology const *> getTopologies() const override;

private:
    struct Impl;
    std::unique_ptr<Impl> pimpl;
    readdy::util::thread::Config const *const config;
};
}
}
}
