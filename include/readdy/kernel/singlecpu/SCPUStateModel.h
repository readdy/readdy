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

#include <readdy/model/KernelStateModel.h>
#include <memory>
#include <readdy/model/Vec3.h>
#include <readdy/kernel/singlecpu/model/SCPUParticleData.h>
#include <readdy/model/KernelContext.h>
#include <readdy/kernel/singlecpu/model/SCPUNeighborList.h>
#include <readdy/model/reactions/ReactionRecord.h>

namespace readdy {
namespace kernel {
namespace scpu {

class SCPUStateModel : public readdy::model::KernelStateModel {
    using topology_action_factory = readdy::model::top::TopologyActionFactory;
public:

    virtual void updateNeighborList() override;

    virtual void clearNeighborList() override;

    virtual void calculateForces() override;

    virtual void addParticle(const readdy::model::Particle &p) override;

    virtual void addParticles(const std::vector<readdy::model::Particle> &p) override;

    virtual void removeParticle(const readdy::model::Particle &p) override;

    virtual void removeAllParticles() override;

    virtual readdy::model::top::GraphTopology *const addTopology(const std::vector<readdy::model::TopologyParticle> &particles) override;

    virtual const std::vector<readdy::model::Vec3> getParticlePositions() const override;

    virtual readdy::model::Particle getParticleForIndex(const std::size_t index) const override;

    virtual double getEnergy() const override;

    virtual void increaseEnergy(double increase);

    SCPUStateModel(readdy::model::KernelContext const *context, const topology_action_factory *const );

    ~SCPUStateModel();

    // move
    SCPUStateModel(SCPUStateModel &&rhs);

    SCPUStateModel &operator=(SCPUStateModel &&rhs);

    virtual readdy::kernel::scpu::model::SCPUParticleData *getParticleData() const;

    virtual const model::SCPUNeighborList *getNeighborList() const;

    virtual const std::vector<readdy::model::Particle> getParticles() const override;

    std::vector<readdy::model::reactions::ReactionRecord>& reactionRecords();

    const std::vector<readdy::model::reactions::ReactionRecord>& reactionRecords() const;

    std::tuple<std::vector<std::size_t>, std::vector<std::size_t>>& reactionCounts();

    const std::tuple<std::vector<std::size_t>, std::vector<std::size_t>>& reactionCounts() const;

    virtual void expected_n_particles(const std::size_t n) override;

private:
    struct Impl;
    std::unique_ptr<Impl> pimpl;
};

}
}
}
