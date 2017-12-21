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
 * @file CPUActionFactory.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 23.06.16
 */

#include <readdy/kernel/cpu_legacy/actions/CPUActionFactory.h>
#include <readdy/kernel/cpu_legacy/actions/CPUEulerBDIntegrator.h>
#include <readdy/kernel/cpu_legacy/actions/CPUUpdateNeighborList.h>
#include <readdy/kernel/cpu_legacy/actions/CPUCalculateForces.h>
#include <readdy/kernel/cpu_legacy/actions/CPUEvaluateCompartments.h>
#include <readdy/kernel/cpu_legacy/actions/reactions/CPUGillespie.h>
#include <readdy/kernel/cpu_legacy/actions/reactions/CPUUncontrolledApproximation.h>
#include <readdy/kernel/cpu_legacy/actions/CPUEvaluateTopologyReactions.h>

namespace core_p = readdy::model::actions;

namespace readdy {
namespace kernel {
namespace cpu_legacy {
namespace actions {
CPUActionFactory::CPUActionFactory(CPULegacyKernel *const kernel) : kernel(kernel) { }

core_p::EulerBDIntegrator *CPUActionFactory::createEulerBDIntegrator(scalar timeStep) const {
    return new CPUEulerBDIntegrator(kernel, timeStep);
}

core_p::CalculateForces *CPUActionFactory::createCalculateForces() const {
    return new CPUCalculateForces(kernel);
}

core_p::UpdateNeighborList *
CPUActionFactory::createUpdateNeighborList(core_p::UpdateNeighborList::Operation operation,
                                            scalar skinSize) const {
    return new CPUUpdateNeighborList(kernel, operation, skinSize);
}

core_p::EvaluateCompartments *CPUActionFactory::createEvaluateCompartments() const {
    return new CPUEvaluateCompartments(kernel);
}

core_p::reactions::UncontrolledApproximation *
CPUActionFactory::createUncontrolledApproximation(scalar timeStep) const {
    return new reactions::CPUUncontrolledApproximation(kernel, timeStep);
}

core_p::reactions::Gillespie *CPUActionFactory::createGillespie(scalar timeStep) const {
    return new reactions::CPUGillespie(kernel, timeStep);
}

readdy::model::actions::AddParticles *
CPUActionFactory::createAddParticles(const std::vector<readdy::model::Particle> &particles) const {
    return new readdy::model::actions::AddParticles(kernel, particles);
}

readdy::model::actions::top::EvaluateTopologyReactions *
CPUActionFactory::createEvaluateTopologyReactions(scalar timeStep) const {
    return new top::CPUEvaluateTopologyReactions(kernel, timeStep);
}
}
}
}
}