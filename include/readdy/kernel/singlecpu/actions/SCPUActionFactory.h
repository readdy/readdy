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


//
// Created by clonker on 08.04.16.
//

#pragma once

#include <readdy/model/actions/ActionFactory.h>
#include <readdy/kernel/singlecpu/SCPUStateModel.h>

namespace readdy {
namespace kernel {
namespace scpu {
class SCPUKernel;
namespace actions {
class SCPUActionFactory : public readdy::model::actions::ActionFactory {
public:
    SCPUActionFactory(SCPUKernel *const kernel);

protected:
    readdy::model::actions::AddParticles *createAddParticles(const std::vector<readdy::model::Particle> &particles) const override;

    readdy::model::actions::EulerBDIntegrator *createEulerBDIntegrator(double timeStep) const override;

    readdy::model::actions::CalculateForces *createCalculateForces() const override;

    readdy::model::actions::UpdateNeighborList *
    createUpdateNeighborList(readdy::model::actions::UpdateNeighborList::Operation, double) const override;

    readdy::model::actions::EvaluateCompartments *createEvaluateCompartments() const override;

    readdy::model::actions::reactions::UncontrolledApproximation *
    createUncontrolledApproximation(double timeStep) const override;

    readdy::model::actions::reactions::Gillespie *createGillespie(double timeStep) const override;

    readdy::model::actions::reactions::GillespieParallel *createGillespieParallel(double timeStep) const override;

    readdy::model::actions::reactions::NextSubvolumes *createNextSubvolumes(double timeStep) const override;

private:
    SCPUKernel *const kernel;
};
}
}
}
}
