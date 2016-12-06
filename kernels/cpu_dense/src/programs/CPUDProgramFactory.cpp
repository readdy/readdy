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
 * @file ProgramFactory.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 23.11.16
 */


#include <readdy/kernel/cpu_dense/programs/CPUDProgramFactory.h>
#include <readdy/kernel/cpu_dense/programs/CPUDEulerBDIntegrator.h>
#include <readdy/kernel/cpu_dense/programs/CPUDUpdateNeighborList.h>
#include <readdy/kernel/cpu_dense/programs/CPUDCalculateForces.h>
#include <readdy/kernel/cpu_dense/programs/CPUDCompartments.h>
#include <readdy/kernel/cpu_dense/programs/reactions/CPUDGillespie.h>
#include <readdy/kernel/cpu_dense/programs/reactions/CPUDUncontrolledApproximation.h>
#include <readdy/kernel/cpu_dense/programs/reactions/CPUDGillespieParallel.h>

namespace core_p = readdy::model::programs;

namespace readdy {
namespace kernel {
namespace cpu_dense {
namespace programs {
CPUDProgramFactory::CPUDProgramFactory(CPUDKernel *kernel) {
    factory[core_p::getProgramName<reactions::CPUDUncontrolledApproximation>()] = [kernel] {
        return new reactions::CPUDUncontrolledApproximation(kernel);
    };
    factory[core_p::getProgramName<CPUDEulerBDIntegrator>()] = [kernel] {
        return new CPUDEulerBDIntegrator(kernel);
    };
    factory[core_p::getProgramName<CPUDUpdateNeighborList>()] = [kernel] {
        return new CPUDUpdateNeighborList(kernel);
    };
    factory[core_p::getProgramName<CPUDCalculateForces>()] = [kernel] {
        return new CPUDCalculateForces(kernel);
    };
    factory[core_p::getProgramName<reactions::CPUDGillespie>()] = [kernel] {
        return new reactions::CPUDGillespie(kernel);
    };
    factory[core_p::getProgramName<reactions::CPUDGillespieParallel>()] = [kernel] {
        return new reactions::CPUDGillespieParallel(kernel);
    };
    factory[core_p::getProgramName<CPUDCompartments>()] = [kernel] {
        return new CPUDCompartments(kernel);
    };
}
}
}
}
}