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
 * @file CPUProgramFactory.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 23.06.16
 */

#include <readdy/kernel/cpu/programs/CPUProgramFactory.h>
#include <readdy/kernel/cpu/programs/CPUEulerBDIntegrator.h>
#include <readdy/kernel/cpu/programs/CPUUpdateNeighborList.h>
#include <readdy/kernel/cpu/programs/CPUCalculateForces.h>
#include <readdy/kernel/cpu/programs/CPUCompartments.h>
#include <readdy/kernel/cpu/programs/reactions/CPUGillespie.h>
#include <readdy/kernel/cpu/programs/reactions/CPUUncontrolledApproximation.h>
#include <readdy/kernel/cpu/programs/reactions/CPUGillespieParallel.h>
#include <readdy/kernel/cpu/programs/reactions/NextSubvolumesReactionScheduler.h>

namespace core_p = readdy::model::actions;

namespace readdy {
namespace kernel {
namespace cpu {
namespace programs {
CPUProgramFactory::CPUProgramFactory(CPUKernel *const kernel) : kernel(kernel) { }

core_p::EulerBDIntegrator *CPUProgramFactory::createEulerBDIntegrator(double timeStep) const {
    return new CPUEulerBDIntegrator(kernel, timeStep);
}

core_p::CalculateForces *CPUProgramFactory::createCalculateForces() const {
    return new CPUCalculateForces(kernel);
}

core_p::UpdateNeighborList *
CPUProgramFactory::createUpdateNeighborList(core_p::UpdateNeighborList::Operation operation,
                                            double skinSize) const {
    return new CPUUpdateNeighborList(kernel, operation, skinSize);
}

core_p::Compartments *CPUProgramFactory::createCompartments() const {
    return new CPUCompartments(kernel);
}

core_p::reactions::UncontrolledApproximation *
CPUProgramFactory::createUncontrolledApproximation(double timeStep) const {
    return new reactions::CPUUncontrolledApproximation(kernel, timeStep);
}

core_p::reactions::Gillespie *CPUProgramFactory::createGillespie(double timeStep) const {
    return new reactions::CPUGillespie(kernel, timeStep);
}

core_p::reactions::GillespieParallel *CPUProgramFactory::createGillespieParallel(double timeStep) const {
    return new reactions::CPUGillespieParallel(kernel, timeStep);
}

core_p::reactions::NextSubvolumes *CPUProgramFactory::createNextSubvolumes(double timeStep) const {
    return new reactions::CPUNextSubvolumes(kernel, timeStep);
}
}
}
}
}