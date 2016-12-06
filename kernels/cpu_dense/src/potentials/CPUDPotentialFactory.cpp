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
 * @file PotentialFactory.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 23.11.16
 */

#include <readdy/kernel/cpu_dense/potentials/CPUDPotentialFactory.h>
#include <readdy/kernel/cpu_dense/CPUDKernel.h>
#include <readdy/kernel/singlecpu/potentials/SCPUPotentialsOrder1.h>
#include <readdy/kernel/singlecpu/potentials/SCPUPotentialsOrder2.h>

namespace readdy {
namespace kernel {
namespace cpu_dense {
namespace potentials {
CPUDPotentialFactory::CPUDPotentialFactory(CPUDKernel *const kernel) {
    namespace p = readdy::model::potentials;
    namespace singlecpu_pot = readdy::kernel::scpu::potentials;
    factory[p::getPotentialName<p::CubePotential>()] = [kernel] {
        return new singlecpu_pot::SCPUCubePotential(kernel);
    };
    factory[p::getPotentialName<p::HarmonicRepulsion>()] = [kernel] {
        return new singlecpu_pot::SCPUHarmonicRepulsion(kernel);
    };
    factory[p::getPotentialName<p::WeakInteractionPiecewiseHarmonic>()] = [kernel] {
        return new singlecpu_pot::SCPUWeakInteractionPiecewiseHarmonic(kernel);
    };
    factory[p::getPotentialName<p::SpherePotential>()] = [kernel] {
        return new singlecpu_pot::SCPUSpherePotential(kernel);
    };
}
}
}
}
}