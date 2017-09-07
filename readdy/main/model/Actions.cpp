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
 * @file Programs.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 17.11.16
 */

#include <readdy/model/actions/Actions.h>
#include <readdy/model/Kernel.h>

namespace readdy {
namespace model {
namespace actions {


UpdateNeighborList::UpdateNeighborList(UpdateNeighborList::Operation operation, scalar skinSize)
        : operation(operation), skinSize(skinSize) {
}

EulerBDIntegrator::EulerBDIntegrator(scalar timeStep) : TimeStepDependentAction(timeStep) {}

reactions::UncontrolledApproximation::UncontrolledApproximation(scalar timeStep) : TimeStepDependentAction(timeStep) {}

reactions::Gillespie::Gillespie(scalar timeStep) : TimeStepDependentAction(timeStep) {}

reactions::GillespieParallel::GillespieParallel(scalar timeStep) : TimeStepDependentAction(timeStep) {}

reactions::NextSubvolumes::NextSubvolumes(scalar timeStep) : TimeStepDependentAction(timeStep) {}

AddParticles::AddParticles(Kernel *const kernel, const std::vector<Particle> &particles)
        : particles(particles), kernel(kernel) {}

AddParticles::AddParticles(Kernel *const kernel, const Particle &particle)
        : AddParticles(kernel, std::vector<Particle>{particle}) {}

void AddParticles::perform(const util::PerformanceNode &node) {
    auto t = node.timeit();
    if(kernel != nullptr) {
        kernel->getKernelStateModel().addParticles(particles);
    } else {
        log::critical("Tried to perform {} without providing a valid kernel!", getActionName<AddParticles>());
    }
}

CalculateForces::CalculateForces() : Action() {}

top::EvaluateTopologyReactions::EvaluateTopologyReactions(scalar timeStep) : TimeStepDependentAction(timeStep) {}

}
}
}