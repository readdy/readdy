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

namespace core_p = readdy::model::programs;

namespace readdy {
namespace kernel {
namespace cpu {
namespace programs {
CPUProgramFactory::CPUProgramFactory(CPUKernel *kernel) {
    factory[core_p::getProgramName<core_p::reactions::UncontrolledApproximation>()] = [kernel] {
        return new reactions::CPUUncontrolledApproximation(kernel);
    };
    factory[core_p::getProgramName<core_p::EulerBDIntegrator>()] = [kernel] {
        return new CPUEulerBDIntegrator(kernel);
    };
    factory[core_p::getProgramName<core_p::UpdateNeighborList>()] = [kernel] {
        return new CPUUpdateNeighborList(kernel);
    };
    factory[core_p::getProgramName<core_p::CalculateForces>()] = [kernel] {
        return new CPUCalculateForces(kernel);
    };
    factory[core_p::getProgramName<core_p::reactions::Gillespie>()] = [kernel] {
        return new reactions::CPUGillespie(kernel);
    };
    factory[core_p::getProgramName<core_p::reactions::GillespieParallel>()] = [kernel] {
        return new reactions::CPUGillespieParallel(kernel);
    };
    factory[core_p::getProgramName<core_p::reactions::NextSubvolumes>()] = [kernel] {
        return new reactions::NextSubvolumes(kernel);
    };
    factory[core_p::getProgramName<core_p::Compartments>()] = [kernel] {
        return new CPUCompartments(kernel);
    };
}
}
}
}
}