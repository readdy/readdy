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
 * @file CPUProgramFactory.h.h
 * @brief << brief description >>
 * @author clonker
 * @date 22.11.16
 */


#ifndef READDY_DENSE_CPUPROGRAMFACTORY_H
#define READDY_DENSE_CPUPROGRAMFACTORY_H

#include <readdy/kernel/cpu_dense/CPUDKernel.h>

namespace readdy {
namespace kernel {
namespace cpu_dense {
namespace programs {
class CPUDProgramFactory : public readdy::model::actions::ActionFactory {
    CPUDKernel *const kernel;
public:
    CPUDProgramFactory(CPUDKernel *const kernel);

    readdy::model::actions::EulerBDIntegrator *createEulerBDIntegrator(double timeStep) const override;

    readdy::model::actions::CalculateForces *createCalculateForces() const override;

    readdy::model::actions::UpdateNeighborList *
    createUpdateNeighborList(readdy::model::actions::UpdateNeighborList::Operation, double) const override;

    readdy::model::actions::Compartments *createCompartments() const override;

    readdy::model::actions::reactions::UncontrolledApproximation *
    createUncontrolledApproximation(double timeStep) const override;

    readdy::model::actions::reactions::Gillespie *createGillespie(double timeStep) const override;

    readdy::model::actions::reactions::GillespieParallel *createGillespieParallel(double timeStep) const override;

    readdy::model::actions::reactions::NextSubvolumes *createNextSubvolumes(double timeStep) const override;
};

}
}
}
}

#endif //READDY_DENSE_CPUPROGRAMFACTORY_H
