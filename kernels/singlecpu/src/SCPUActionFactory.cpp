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

#include <memory>

#include <readdy/kernel/singlecpu/actions/SCPUActionFactory.h>
#include <readdy/kernel/singlecpu/actions/SCPUEulerBDIntegrator.h>
#include <readdy/kernel/singlecpu/actions/SCPUCalculateForces.h>
#include <readdy/kernel/singlecpu/actions/SCPUReactionImpls.h>
#include <readdy/kernel/singlecpu/actions/SCPUUpdateNeighborList.h>
#include <readdy/kernel/singlecpu/actions/SCPUEvaluateCompartments.h>
#include <readdy/kernel/singlecpu/actions/SCPUEvaluateTopologyReactions.h>

namespace core_actions = readdy::model::actions;

namespace readdy {
namespace kernel {
namespace scpu {
namespace actions {
SCPUActionFactory::SCPUActionFactory(SCPUKernel *const kernel) : kernel(kernel) {}

namespace rma = readdy::model::actions;

std::vector<std::string> SCPUActionFactory::getAvailableActions() const {
    return {
            rma::getActionName<rma::AddParticles>(), rma::getActionName<rma::EulerBDIntegrator>(),
            rma::getActionName<rma::CalculateForces>(),
            rma::getActionName<rma::UpdateNeighborList>(),
            rma::getActionName<rma::reactions::UncontrolledApproximation>(),
            rma::getActionName<rma::reactions::Gillespie>(),
            rma::getActionName<rma::top::EvaluateTopologyReactions>()
    };
}

std::unique_ptr<readdy::model::actions::EulerBDIntegrator> SCPUActionFactory::eulerBDIntegrator(scalar timeStep) const {
    return {std::make_unique<SCPUEulerBDIntegrator>(kernel, timeStep)};
}

std::unique_ptr<readdy::model::actions::CalculateForces> SCPUActionFactory::calculateForces() const {
    return {std::make_unique<SCPUCalculateForces>(kernel)};
}

std::unique_ptr<readdy::model::actions::UpdateNeighborList>
SCPUActionFactory::updateNeighborList(readdy::model::actions::UpdateNeighborList::Operation operation,
                                      scalar skinSize) const {
    return {std::make_unique<SCPUUpdateNeighborList>(kernel, operation, skinSize)};
}

std::unique_ptr<readdy::model::actions::EvaluateCompartments> SCPUActionFactory::evaluateCompartments() const {
    return {std::make_unique<SCPUEvaluateCompartments>(kernel)};
}

std::unique_ptr<readdy::model::actions::reactions::UncontrolledApproximation>
SCPUActionFactory::uncontrolledApproximation(scalar timeStep) const {
    return {std::make_unique<reactions::SCPUUncontrolledApproximation>(kernel, timeStep)};
}

std::unique_ptr<readdy::model::actions::reactions::Gillespie> SCPUActionFactory::gillespie(scalar timeStep) const {
    return {std::make_unique<reactions::SCPUGillespie>(kernel, timeStep)};
}

std::unique_ptr<readdy::model::actions::top::EvaluateTopologyReactions>
SCPUActionFactory::evaluateTopologyReactions(scalar timeStep) const {
    return {std::make_unique<top::SCPUEvaluateTopologyReactions>(kernel, timeStep)};
}

std::unique_ptr<readdy::model::actions::reactions::DetailedBalance>
SCPUActionFactory::detailedBalance(scalar timeStep) const {
    return {std::make_unique<reactions::SCPUDetailedBalance>(kernel, timeStep)};
}

}
}
}
}


