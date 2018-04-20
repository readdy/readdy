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


#include <readdy/kernel/singlecpu/SCPUKernel.h>

#include <readdy/kernel/singlecpu/observables/SCPUObservableFactory.h>
#include <readdy/kernel/singlecpu/observables/SCPUObservables.h>
#include <readdy/kernel/singlecpu/observables/SCPUAggregators.h>

namespace readdy {
namespace kernel {
namespace scpu {
namespace observables {
SCPUObservableFactory::SCPUObservableFactory(readdy::kernel::scpu::SCPUKernel *const kernel)
        : ObservableFactory(kernel), kernel(kernel) {
}

std::unique_ptr<readdy::model::observables::HistogramAlongAxis>
SCPUObservableFactory::histogramAlongAxis(stride_type stride, std::vector<scalar> binBorders,
                                          std::vector<std::string> typesToCount, unsigned int axis) const {
    return {std::make_unique<SCPUHistogramAlongAxis>(kernel, stride, binBorders, typesToCount, axis)};
}

std::unique_ptr<readdy::model::observables::NParticles>
SCPUObservableFactory::nParticles(stride_type stride, std::vector<std::string> typesToCount) const {
    return {std::make_unique<SCPUNParticles>(kernel, stride, typesToCount)};
}

std::unique_ptr<readdy::model::observables::Forces>
SCPUObservableFactory::forces(stride_type stride, std::vector<std::string> typesToCount) const {
    return {std::make_unique<SCPUForces>(kernel, stride, typesToCount)};
}

std::unique_ptr<readdy::model::observables::Positions>
SCPUObservableFactory::positions(stride_type stride, std::vector<std::string> typesToCount) const {
    return {std::make_unique<SCPUPositions>(kernel, stride, typesToCount)};
}

std::unique_ptr<readdy::model::observables::RadialDistribution>
SCPUObservableFactory::radialDistribution(stride_type stride, std::vector<scalar> binBorders,
                                          std::vector<std::string> typeCountFrom, std::vector<std::string> typeCountTo,
                                          scalar particleDensity) const {
    return {std::make_unique<readdy::model::observables::RadialDistribution>(kernel, stride, binBorders, typeCountFrom,
                                                                             typeCountTo, particleDensity)};
}

std::unique_ptr<readdy::model::observables::Particles>
SCPUObservableFactory::particles(stride_type stride) const {
    return {std::make_unique<SCPUParticles>(kernel, stride)};
}

std::unique_ptr<readdy::model::observables::MeanSquaredDisplacement>
SCPUObservableFactory::msd(stride_type stride, std::vector<std::string> typesToCount,
                           readdy::model::observables::Particles *particlesObservable) const {
    return {std::make_unique<SCPUMeanSquaredDisplacement<>>(kernel, stride, typesToCount, particlesObservable)};
}

std::unique_ptr<readdy::model::observables::Reactions>
SCPUObservableFactory::reactions(readdy::model::observables::ObservableFactory::stride_type stride) const {
    return {std::make_unique<SCPUReactions>(kernel, stride)};
}

std::unique_ptr<readdy::model::observables::ReactionCounts>
SCPUObservableFactory::reactionCounts(readdy::model::observables::ObservableFactory::stride_type stride) const {
    return {std::make_unique<SCPUReactionCounts>(kernel, stride)};
}

std::unique_ptr<readdy::model::observables::Virial>
SCPUObservableFactory::virial(readdy::model::observables::ObservableFactory::stride_type stride) const {
    return {std::make_unique<SCPUVirial>(kernel, stride)};
}

}
}
}
}