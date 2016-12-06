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
 * @file SingleCPUEulerDBIntegrator.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 19.04.16
 */

#include <readdy/kernel/singlecpu/programs/SCPUEulerBDIntegrator.h>

namespace readdy {
namespace kernel {
namespace scpu {
namespace programs {

SCPUEulerBDIntegrator::SCPUEulerBDIntegrator(SCPUKernel *kernel)
        : readdy::model::programs::EulerBDIntegrator(), kernel(kernel) {};

void SCPUEulerBDIntegrator::execute() {
    const auto &context = kernel->getKernelContext();
    const auto &kbt = context.getKBT();
    const auto &fixPos = context.getFixPositionFun();
    const auto &&dt = context.getTimeStep();
    const auto &&pd = kernel->getKernelStateModel().getParticleData();
    auto it_pos = pd->begin_positions();
    auto it_types = pd->begin_types();
    auto it_forces = pd->begin_forces();
    for (; it_pos != pd->end_positions();) {
        const double D = context.getDiffusionConstant(*it_types);
        const auto randomDisplacement = std::sqrt(2. * D * dt) * (readdy::model::rnd::normal3());
        *it_pos += randomDisplacement;
        const auto deterministicDisplacement = *it_forces * dt * D / kbt;
        *it_pos += deterministicDisplacement;
        fixPos(*it_pos);
        ++it_pos;
        ++it_types;
        ++it_forces;
    }
}
}
}
}
}


