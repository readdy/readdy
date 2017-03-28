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
#include <readdy/kernel/cpu/model/CPUNeighborList.h>
#include <readdy/common/thread/Config.h>
#include <readdy/kernel/cpu/model/CPUParticleData.h>
#include <readdy/model/reactions/ReactionRecord.h>

namespace readdy {
namespace kernel {
namespace cpu {
class CPUStateModel : public readdy::model::KernelStateModel {

public:

    using data_t = readdy::kernel::cpu::model::CPUParticleData;

    CPUStateModel(readdy::model::KernelContext *const context, readdy::util::thread::Config const *const config,
                  readdy::model::top::TopologyActionFactory const *const taf);

    ~CPUStateModel();

    virtual const std::vector<readdy::model::Vec3> getParticlePositions() const override;

    virtual const std::vector<readdy::model::Particle> getParticles() const override;

    virtual void updateNeighborList() override;

    virtual void calculateForces() override;

    virtual void addParticle(const readdy::model::Particle &p) override;

    virtual void addParticles(const std::vector<readdy::model::Particle> &p) override;

    virtual void removeParticle(const readdy::model::Particle &p) override;

    virtual void removeAllParticles() override;

    virtual double getEnergy() const override;

    data_t const *const getParticleData() const;

    data_t *const getParticleData();

    model::CPUNeighborList const *const getNeighborList() const;

    model::CPUNeighborList *const getNeighborList();

    virtual void expected_n_particles(const std::size_t n) override;

    virtual void clearNeighborList() override;

    virtual readdy::model::top::GraphTopology *const
    addTopology(const std::vector<readdy::model::TopologyParticle> &particles) override;

    std::vector<readdy::model::reactions::ReactionRecord> &reactionRecords();

    const std::vector<readdy::model::reactions::ReactionRecord> &reactionRecords() const;

    std::tuple<std::vector<std::size_t>, std::vector<std::size_t>> &reactionCounts();

    const std::tuple<std::vector<std::size_t>, std::vector<std::size_t>> &reactionCounts() const;

    virtual readdy::model::Particle getParticleForIndex(const std::size_t index) const override;

private:
    struct Impl;
    std::unique_ptr<Impl> pimpl;
    readdy::util::thread::Config const *const config;
};
}
}
}
