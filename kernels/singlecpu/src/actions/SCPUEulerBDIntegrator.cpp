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

#include <readdy/kernel/singlecpu/actions/SCPUEulerBDIntegrator.h>

namespace readdy {
namespace kernel {
namespace scpu {
namespace actions {

SCPUEulerBDIntegrator::SCPUEulerBDIntegrator(SCPUKernel *kernel, double timeStep)
        : readdy::model::actions::EulerBDIntegrator(timeStep), kernel(kernel) {};

void SCPUEulerBDIntegrator::perform() {
    const auto &context = kernel->getKernelContext();
    const auto &kbt = context.getKBT();
    const auto &fixPos = context.getFixPositionFun();
    auto& stateModel = kernel->getSCPUKernelStateModel();
    const auto pd = stateModel.getParticleData();
    for(auto& entry : *pd) {
        if(!entry.is_deactivated()) {
            const double D = context.particle_types().diffusion_constant_of(entry.type);
            const auto randomDisplacement = std::sqrt(2. * D * timeStep) * (readdy::model::rnd::normal3());
            entry.pos += randomDisplacement;
            const auto deterministicDisplacement = entry.force * timeStep * D / kbt;
            entry.pos += deterministicDisplacement;
            fixPos(entry.pos);
        }
    }
}
}
}
}
}


