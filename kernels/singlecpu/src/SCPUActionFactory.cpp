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

#include <readdy/common/make_unique.h>
#include <readdy/kernel/singlecpu/programs/SCPUActionFactory.h>
#include <readdy/kernel/singlecpu/programs/SCPUEulerBDIntegrator.h>
#include <readdy/kernel/singlecpu/programs/SCPUCalculateForces.h>
#include <readdy/kernel/singlecpu/programs/SCPUReactionImpls.h>
#include <readdy/kernel/singlecpu/programs/SCPUUpdateNeighborList.h>
#include <readdy/kernel/singlecpu/programs/SCPUCompartments.h>

namespace core_actions = readdy::model::actions;

namespace readdy {
namespace kernel {
namespace scpu {
namespace actions {
SCPUActionFactory::SCPUActionFactory(SCPUKernel *const kernel) : kernel(kernel) {}

core_actions::EulerBDIntegrator *SCPUActionFactory::createEulerBDIntegrator(double timeStep) const {
    return new SCPUEulerBDIntegrator(kernel, timeStep);
}

core_actions::CalculateForces *SCPUActionFactory::createCalculateForces() const {
    return new SCPUCalculateForces(kernel);
}

core_actions::UpdateNeighborList *
SCPUActionFactory::createUpdateNeighborList(core_actions::UpdateNeighborList::Operation op, double skinSize) const {
    return new SCPUUpdateNeighborList(kernel, op, skinSize);
}

core_actions::Compartments *SCPUActionFactory::createCompartments() const {
    return new SCPUCompartments(kernel);
}

core_actions::reactions::UncontrolledApproximation *
SCPUActionFactory::createUncontrolledApproximation(double timeStep) const {
    return new reactions::SCPUUncontrolledApproximation(kernel, timeStep);
}

core_actions::reactions::Gillespie *SCPUActionFactory::createGillespie(double timeStep) const {
    return new reactions::SCPUGillespie(kernel, timeStep);
}

core_actions::reactions::GillespieParallel *SCPUActionFactory::createGillespieParallel(double) const {
    log::console()->critical("SingleCPU kernel does not support the \"{}\" action",
                             core_actions::getActionName<core_actions::reactions::GillespieParallel>());
    return nullptr;
}

core_actions::reactions::NextSubvolumes *SCPUActionFactory::createNextSubvolumes(double) const {
    log::console()->critical("SingleCPU kernel does not support the \"{}\" action",
                             core_actions::getActionName<core_actions::reactions::NextSubvolumes>());
    return nullptr;
}

}
}
}
}


