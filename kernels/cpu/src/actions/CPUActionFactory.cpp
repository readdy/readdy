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

#include <readdy/kernel/cpu/actions/CPUActionFactory.h>
#include <readdy/kernel/cpu/actions/CPUEulerBDIntegrator.h>
#include <readdy/kernel/cpu/actions/CPUUpdateNeighborList.h>
#include <readdy/kernel/cpu/actions/CPUCalculateForces.h>
#include <readdy/kernel/cpu/actions/CPUEvaluateCompartments.h>
#include <readdy/kernel/cpu/actions/reactions/CPUGillespie.h>
#include <readdy/kernel/cpu/actions/reactions/CPUUncontrolledApproximation.h>
#include <readdy/kernel/cpu/actions/CPUEvaluateTopologyReactions.h>

namespace core_p = readdy::model::actions;

namespace readdy {
namespace kernel {
namespace cpu {
namespace actions {
CPUActionFactory::CPUActionFactory(CPUKernel *const kernel) : kernel(kernel) { }

std::unique_ptr<model::actions::AddParticles>
CPUActionFactory::addParticles(const std::vector<model::Particle> &particles) const {
    return {std::make_unique<readdy::model::actions::AddParticles>(kernel, particles)};
}

std::unique_ptr<model::actions::EulerBDIntegrator> CPUActionFactory::eulerBDIntegrator(scalar timeStep) const {
    return {std::make_unique<CPUEulerBDIntegrator>(kernel, timeStep)};
}

std::unique_ptr<model::actions::CalculateForces> CPUActionFactory::calculateForces() const {
    return {std::make_unique<CPUCalculateForces>(kernel)};
}

std::unique_ptr<model::actions::UpdateNeighborList>
CPUActionFactory::updateNeighborList(model::actions::UpdateNeighborList::Operation operation, scalar skinSize) const {
    return {std::make_unique<CPUUpdateNeighborList>(kernel, operation, skinSize)};
}

std::unique_ptr<model::actions::EvaluateCompartments> CPUActionFactory::evaluateCompartments() const {
    return {std::make_unique<CPUEvaluateCompartments>(kernel)};
}

std::unique_ptr<model::actions::reactions::UncontrolledApproximation>
CPUActionFactory::uncontrolledApproximation(scalar timeStep) const {
    return {std::make_unique<reactions::CPUUncontrolledApproximation>(kernel, timeStep)};
}

std::unique_ptr<model::actions::reactions::Gillespie> CPUActionFactory::gillespie(scalar timeStep) const {
    return {std::make_unique<reactions::CPUGillespie>(kernel, timeStep)};
}

std::unique_ptr<model::actions::top::EvaluateTopologyReactions>
CPUActionFactory::evaluateTopologyReactions(scalar timeStep) const {
    return {std::make_unique<top::CPUEvaluateTopologyReactions>(kernel, timeStep)};
}

std::unique_ptr<readdy::model::actions::reactions::DetailedBalance>
CPUActionFactory::detailedBalance(scalar timeStep) const {
    throw std::invalid_argument("DetailedBalance reaction handler not implemented for CPU");
}
}
}
}
}
