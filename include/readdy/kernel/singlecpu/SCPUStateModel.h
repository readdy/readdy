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

#ifndef READDY_MAIN_SINGLECPUKERNELSTATEMODEL_H
#define READDY_MAIN_SINGLECPUKERNELSTATEMODEL_H


#include <readdy/model/KernelStateModel.h>
#include <memory>
#include <readdy/model/Vec3.h>
#include <readdy/kernel/singlecpu/model/SCPUParticleData.h>
#include <readdy/model/KernelContext.h>
#include <readdy/kernel/singlecpu/model/SCPUNeighborList.h>

namespace readdy {
namespace kernel {
namespace scpu {

class SCPUStateModel : public readdy::model::KernelStateModel {
public:

    virtual void updateNeighborList() override;

    virtual void clearNeighborList() override;

    virtual void calculateForces() override;

    virtual void addParticle(const readdy::model::Particle &p) override;

    virtual void addParticles(const std::vector<readdy::model::Particle> &p) override;

    virtual void removeParticle(const readdy::model::Particle &p) override;

    virtual void removeAllParticles() override;

    virtual const std::vector<readdy::model::Vec3> getParticlePositions() const override;

    virtual double getEnergy() const override;

    virtual void increaseEnergy(double increase);

    SCPUStateModel(readdy::model::KernelContext const *context);

    ~SCPUStateModel();

    // move
    SCPUStateModel(SCPUStateModel &&rhs);

    SCPUStateModel &operator=(SCPUStateModel &&rhs);

    virtual readdy::kernel::scpu::model::SCPUParticleData *getParticleData() const;

    virtual const model::SCPUNeighborList *getNeighborList() const;

    virtual const std::vector<readdy::model::Particle> getParticles() const override;

private:
    struct Impl;
    std::unique_ptr<Impl> pimpl;
};

}
}
}

#endif //READDY_MAIN_SINGLECPUKERNELSTATEMODEL_H
