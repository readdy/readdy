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


bool UpdateNeighborList::supportsSkin() const {
    return false;
}

UpdateNeighborList::UpdateNeighborList(UpdateNeighborList::Operation operation, double skinSize)
        : operation(operation), skinSize(skinSize) {
    if(skinSize >= 0 && !supportsSkin()) {
        log::console()->warn("The selected kernel has no Verlet list implemented, thus ignoring the skin size");
    }
}

EulerBDIntegrator::EulerBDIntegrator(double timeStep) : TimeStepDependentAction(timeStep) {}

reactions::UncontrolledApproximation::UncontrolledApproximation(double timeStep) : TimeStepDependentAction(timeStep) {}

reactions::Gillespie::Gillespie(double timeStep) : TimeStepDependentAction(timeStep) {}

reactions::GillespieParallel::GillespieParallel(double timeStep) : TimeStepDependentAction(timeStep) {}

reactions::NextSubvolumes::NextSubvolumes(double timeStep) : TimeStepDependentAction(timeStep) {}

AddParticles::AddParticles(Kernel *const kernel, const std::vector<Particle> &particles)
        : particles(particles), kernel(kernel) {}

AddParticles::AddParticles(Kernel *const kernel, const Particle &particle)
        : AddParticles(kernel, std::vector<Particle>{particle}) {}

void AddParticles::perform() {
    if(kernel) {
        kernel->getKernelStateModel().addParticles(particles);
    } else {
        log::console()->critical("Tried to perform {} without providing a valid kernel!", getActionName<AddParticles>());
    }
}

CalculateForces::CalculateForces() : Action() {}
}
}
}