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

std::unique_ptr<model::observables::HistogramAlongAxis>
CPUObservableFactory::histogramAlongAxis(stride_type stride, std::vector<scalar> binBorders,
                                         std::vector<std::string> typesToCount, unsigned int axis) const {
    return {std::make_unique<CPUHistogramAlongAxis>(kernel, stride, binBorders, typesToCount, axis)};
}

std::unique_ptr<model::observables::NParticles>
CPUObservableFactory::nParticles(stride_type stride, std::vector<std::string> typesToCount) const {
    return {std::make_unique<CPUNParticles>(kernel, stride, typesToCount)};
}

std::unique_ptr<model::observables::Forces>
CPUObservableFactory::forces(stride_type stride, std::vector<std::string> typesToCount) const {
    return {std::make_unique<CPUForces>(kernel, stride, typesToCount)};
}

std::unique_ptr<model::observables::Positions>
CPUObservableFactory::positions(stride_type stride, std::vector<std::string> typesToCount) const {
    return {std::make_unique<CPUPositions>(kernel, stride, typesToCount)};
}

std::unique_ptr<model::observables::RadialDistribution>
CPUObservableFactory::radialDistribution(stride_type stride, std::vector<scalar> binBorders,
                                         std::vector<std::string> typeCountFrom, std::vector<std::string> typeCountTo,
                                         scalar particleDensity) const {
    return {std::make_unique<model::observables::RadialDistribution>(
            kernel, stride, binBorders, typeCountFrom, typeCountTo, particleDensity
    )};
}

std::unique_ptr<model::observables::Particles> CPUObservableFactory::particles(stride_type stride) const {
    return {std::make_unique<CPUParticles>(kernel, stride)};
}

std::unique_ptr<model::observables::MeanSquaredDisplacement>
CPUObservableFactory::msd(stride_type stride, std::vector<std::string> typesToCount,
                          model::observables::Particles *particlesObservable) const {
    return {std::make_unique<readdy::kernel::scpu::observables::SCPUMeanSquaredDisplacement<CPUKernel>>(
            kernel, stride, typesToCount, particlesObservable
    )};
}

std::unique_ptr<model::observables::Reactions> CPUObservableFactory::reactions(stride_type stride) const {
    return {std::make_unique<CPUReactions>(kernel, stride)};
}

std::unique_ptr<model::observables::ReactionCounts> CPUObservableFactory::reactionCounts(stride_type stride) const {
    return {std::make_unique<CPUReactionCounts>(kernel, stride)};
}

std::unique_ptr<model::observables::Virial>
CPUObservableFactory::virial(stride_type stride) const {
    return {std::make_unique<CPUVirial>(kernel, stride)};
}

}
}
}
}
