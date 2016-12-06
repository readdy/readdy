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
 * @file SingleCPUPotentialFactory.cpp
 * @brief Fill the SingleCPUPotentialFactory with kernel-specific constructors for potentials.
 * @author clonker
 * @date 09.06.16
 */

#include <readdy/kernel/singlecpu/SCPUKernel.h>
#include <readdy/kernel/singlecpu/potentials/SCPUPotentialFactory.h>
#include <readdy/kernel/singlecpu/potentials/SCPUPotentialsOrder1.h>
#include <readdy/kernel/singlecpu/potentials/SCPUPotentialsOrder2.h>

namespace readdy {
namespace kernel {
namespace scpu {
namespace potentials {
SCPUPotentialFactory::SCPUPotentialFactory(SCPUKernel *const kernel)
        : readdy::model::potentials::PotentialFactory(), kernel(kernel) {
    namespace p = readdy::model::potentials;
    // order 1
    factory[p::getPotentialName<p::CubePotential>()] = [kernel] { return new SCPUCubePotential(kernel); };
    factory[p::getPotentialName<p::SpherePotential>()] = [kernel] { return new SCPUSpherePotential(kernel); };
    // order 2
    factory[p::getPotentialName<p::HarmonicRepulsion>()] = [kernel] { return new SCPUHarmonicRepulsion(kernel); };
    factory[p::getPotentialName<p::WeakInteractionPiecewiseHarmonic>()] = [kernel] {
        return new SCPUWeakInteractionPiecewiseHarmonic(kernel);
    };
}

}
}
}
}