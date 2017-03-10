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
 * @date 23.11.16
 */


#include <readdy/kernel/cpu_dense/observables/CPUDObservableFactory.h>
#include <readdy/kernel/cpu_dense/CPUDKernel.h>
#include <readdy/kernel/cpu_dense/observables/CPUDObservables.h>
#include <readdy/kernel/singlecpu/observables/SCPUObservables.h>
#include <readdy/kernel/singlecpu/observables/SCPUAggregators.h>

namespace readdy {
namespace kernel {
namespace cpu_dense {
namespace observables {
CPUDObservableFactory::CPUDObservableFactory(CPUDKernel *const kernel) : readdy::model::observables::ObservableFactory(kernel),
                                                             kernel(kernel) {
}

readdy::model::observables::NParticles *
CPUDObservableFactory::createNParticles(unsigned int stride, std::vector<std::string> typesToCount) const {
    return new CPUDNParticles(kernel, stride, typesToCount);
}

readdy::model::observables::HistogramAlongAxis *
CPUDObservableFactory::createHistogramAlongAxis(unsigned int stride, std::vector<double> binBorders,
                                                std::vector<std::string> typesToCount,
                                                unsigned int axis) const {
    return new CPUDHistogramAlongAxis(kernel, stride, binBorders, typesToCount, axis);
}

readdy::model::observables::Forces *
CPUDObservableFactory::createForces(unsigned int stride, std::vector<std::string> typesToCount) const {
    return new CPUDForces(kernel, stride, typesToCount);
}

readdy::model::observables::Positions *
CPUDObservableFactory::createPositions(unsigned int stride, std::vector<std::string> typesToCount) const {
    return new CPUDPositions(kernel, stride, typesToCount);
}

readdy::model::observables::Particles *CPUDObservableFactory::createParticles(unsigned int stride) const {
    return new CPUDParticles(kernel, stride);
}

readdy::model::observables::RadialDistribution *
CPUDObservableFactory::createRadialDistribution(unsigned int stride, std::vector<double> binBorders, std::vector<std::string> typeCountFrom,
                                                std::vector<std::string> typeCountTo, double particleToDensity) const {
    return new readdy::kernel::scpu::observables::SCPURadialDistribution<CPUDKernel>(kernel, stride, binBorders, typeCountFrom, typeCountTo,
                                                                                    particleToDensity);
}

readdy::model::observables::MeanSquaredDisplacement *
CPUDObservableFactory::createMeanSquaredDisplacement(unsigned int stride, std::vector<std::string> typesToCount,
                                                     readdy::model::observables::Particles *particlesObservable) const {
    return new readdy::kernel::scpu::observables::SCPUMeanSquaredDisplacement<CPUDKernel>(kernel, stride, typesToCount, particlesObservable);
}

readdy::model::observables::Reactions *CPUDObservableFactory::createReactions(unsigned int stride, bool recordPosition) const {
    throw std::runtime_error("Reactions observable is not supported on the CPU_Dense kernel.");
}

}
}
}
}