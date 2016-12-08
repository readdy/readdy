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
 * @file SingleCPUObservableFactory.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 30.06.16
 */


#include <readdy/model/Kernel.h>
#include <readdy/kernel/singlecpu/observables/SCPUObservableFactory.h>
#include <readdy/kernel/singlecpu/observables/SCPUObservables.h>

namespace readdy {
namespace kernel {
namespace scpu {
namespace observables {
SCPUObservableFactory::SCPUObservableFactory(readdy::kernel::scpu::SCPUKernel *const kernel)
        : ObservableFactory(kernel), kernel(kernel) {
}

readdy::model::observables::HistogramAlongAxis *
SCPUObservableFactory::createHistogramAlongAxis(unsigned int stride, std::vector<double> binBorders,
                                                std::vector<std::string> typesToCount,
                                                unsigned int axis) const {
    return new SCPUHistogramAlongAxis(kernel, stride, binBorders, typesToCount, axis);
}

readdy::model::observables::NParticles *SCPUObservableFactory::createNParticles(
        unsigned int stride, std::vector<std::string> typesToCount) const {
    return new SCPUNParticles(kernel, stride, typesToCount);
}

readdy::model::observables::Forces *
SCPUObservableFactory::createForces(unsigned int stride, std::vector<std::string> typesToCount) const {
    return new SCPUForces(kernel, stride, typesToCount);
}

readdy::model::observables::Positions *
SCPUObservableFactory::createPositions(unsigned int stride, std::vector<std::string> typesToCount) const {
    return new SCPUPositions(kernel, stride, typesToCount);
}

readdy::model::observables::RadialDistribution *
SCPUObservableFactory::createRadialDistribution(unsigned int stride, std::vector<double> binBorders, std::vector<std::string> typeCountFrom,
                                                std::vector<std::string> typeCountTo, double particleToDensity) const {
    return new SCPURadialDistribution<SCPUKernel>(kernel, stride, binBorders, typeCountFrom, typeCountTo, particleToDensity);
}

readdy::model::observables::Particles *SCPUObservableFactory::createParticles(unsigned int stride) const {
    return new SCPUParticles(kernel, stride);
}
}
}
}
}