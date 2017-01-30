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
 * @file SCPUTopologyActions.h
 * @brief << brief description >>
 * @author clonker
 * @date 30.01.17
 * @copyright GNU Lesser General Public License v3.0
 */

#ifndef READDY_MAIN_SCPUTOPOLOGYACTIONS_H
#define READDY_MAIN_SCPUTOPOLOGYACTIONS_H

#include <readdy/common/macros.h>
#include <readdy/model/topologies/actions/TopologyActions.h>
#include <readdy/kernel/singlecpu/SCPUStateModel.h>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(kernel)
NAMESPACE_BEGIN(scpu)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(top)

class SCPUCalculateHarmonicBondPotential : public readdy::model::top::CalculateHarmonicBondPotential {

    const readdy::model::top::HarmonicBondPotential *const potential;
    const SCPUParticleData *const data;

public:
    SCPUCalculateHarmonicBondPotential(const readdy::model::KernelContext *const context,
                                       const SCPUParticleData *const data,
                                       const readdy::model::top::HarmonicBondPotential *const potential)
            : CalculateHarmonicBondPotential(context), potential(potential), data(data) {}

    virtual double calculateForcesAndEnergy() override {
        readdy::model::Vec3::entry_t energy = 0;
        const auto& d = context->getShortestDifferenceFun();
        auto itBeginPos = data->begin_positions();
        std::vector<readdy::model::Vec3>::iterator itBeginForce = data->begin_forces();
        for(const auto& bond : potential->getBonds()) {
            const auto x_ij = d(*(itBeginPos + bond.idx1), *(itBeginPos + bond.idx2));
            potential->calculateForce(*(itBeginForce + bond.idx1), x_ij, bond);
            energy += potential->calculateEnergy(x_ij, bond);
        }
        return energy;
    }

};

NAMESPACE_END(top)
NAMESPACE_END(model)
NAMESPACE_END(scpu)
NAMESPACE_END(kernel)
NAMESPACE_END(readdy)
#endif //READDY_MAIN_SCPUTOPOLOGYACTIONS_H
