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
 * @file CPUProgramFactory.h
 * @brief << brief description >>
 * @author clonker
 * @date 23.06.16
 */

#pragma once
#include <readdy/kernel/cpu/CPUKernel.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace actions {
class CPUActionFactory : public readdy::model::actions::ActionFactory {
    CPUKernel *const kernel;
public:
    explicit CPUActionFactory(CPUKernel* kernel);

protected:
    readdy::model::actions::AddParticles *createAddParticles(const std::vector<readdy::model::Particle> &particles) const override;

    readdy::model::actions::EulerBDIntegrator *createEulerBDIntegrator(readdy::scalar timeStep) const override;

    readdy::model::actions::CalculateForces *createCalculateForces() const override;

    readdy::model::actions::UpdateNeighborList *
    createUpdateNeighborList(readdy::model::actions::UpdateNeighborList::Operation operation, readdy::scalar skinSize) const override;

    readdy::model::actions::EvaluateCompartments *createEvaluateCompartments() const override;

    readdy::model::actions::reactions::UncontrolledApproximation *
    createUncontrolledApproximation(readdy::scalar timeStep) const override;

    readdy::model::actions::reactions::Gillespie *createGillespie(readdy::scalar timeStep) const override;

    readdy::model::actions::top::EvaluateTopologyReactions *
    createEvaluateTopologyReactions(readdy::scalar timeStep) const override;
};

}
}
}
}
