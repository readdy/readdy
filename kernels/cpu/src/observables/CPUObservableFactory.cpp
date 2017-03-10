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
 * @file ObservableFactory.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 21.07.16
 */

#include <readdy/kernel/cpu/observables/CPUObservableFactory.h>
#include <readdy/kernel/cpu/CPUKernel.h>
#include <readdy/kernel/cpu/observables/CPUObservables.h>
#include <readdy/kernel/singlecpu/observables/SCPUObservables.h>
#include <readdy/kernel/singlecpu/observables/SCPUAggregators.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace observables {
CPUObservableFactory::CPUObservableFactory(CPUKernel *const kernel) : readdy::model::observables::ObservableFactory(kernel),
                                                                kernel(kernel) {
}

readdy::model::observables::NParticles *
CPUObservableFactory::createNParticles(unsigned int stride, std::vector<std::string> typesToCount) const {
    return new CPUNParticles(kernel, stride, typesToCount);
}

readdy::model::observables::HistogramAlongAxis *
CPUObservableFactory::createHistogramAlongAxis(unsigned int stride, std::vector<double> binBorders,
                                               std::vector<std::string> typesToCount,
                                               unsigned int axis) const {
    return new CPUHistogramAlongAxis(kernel, stride, binBorders, typesToCount, axis);
}

readdy::model::observables::Forces *
CPUObservableFactory::createForces(unsigned int stride, std::vector<std::string> typesToCount) const {
    return new CPUForces(kernel, stride, typesToCount);
}

readdy::model::observables::Positions *
CPUObservableFactory::createPositions(unsigned int stride, std::vector<std::string> typesToCount) const {
    return new CPUPositions(kernel, stride, typesToCount);
}

readdy::model::observables::RadialDistribution *
CPUObservableFactory::createRadialDistribution(unsigned int stride, std::vector<double> binBorders, std::vector<std::string> typeCountFrom,
                                               std::vector<std::string> typeCountTo, double particleToDensity) const {
    return new readdy::kernel::scpu::observables::SCPURadialDistribution<CPUKernel>(kernel, stride, binBorders, typeCountFrom, typeCountTo,
                                                                                               particleToDensity);
}

readdy::model::observables::Particles *CPUObservableFactory::createParticles(unsigned int stride) const {
    return new CPUParticles(kernel, stride);
}

readdy::model::observables::MeanSquaredDisplacement *
CPUObservableFactory::createMeanSquaredDisplacement(unsigned int stride, std::vector<std::string> typesToCount,
                                                    readdy::model::observables::Particles *particlesObservable) const {
    return new readdy::kernel::scpu::observables::SCPUMeanSquaredDisplacement<CPUKernel>(kernel, stride, typesToCount, particlesObservable);
}

readdy::model::observables::Reactions *CPUObservableFactory::createReactions(unsigned int stride, bool recordPosition) const {
    return new CPUReactions(kernel, stride, recordPosition);
}

}
}
}
}
