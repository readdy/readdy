/********************************************************************
 * Copyright © 2018 Computational Molecular Biology Group,          *
 *                  Freie Universität Berlin (GER)                  *
 *                                                                  *
 * Redistribution and use in source and binary forms, with or       *
 * without modification, are permitted provided that the            *
 * following conditions are met:                                    *
 *  1. Redistributions of source code must retain the above         *
 *     copyright notice, this list of conditions and the            *
 *     following disclaimer.                                        *
 *  2. Redistributions in binary form must reproduce the above      *
 *     copyright notice, this list of conditions and the following  *
 *     disclaimer in the documentation and/or other materials       *
 *     provided with the distribution.                              *
 *  3. Neither the name of the copyright holder nor the names of    *
 *     its contributors may be used to endorse or promote products  *
 *     derived from this software without specific                  *
 *     prior written permission.                                    *
 *                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND           *
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,      *
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF         *
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE         *
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR            *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,     *
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,         *
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER *
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,      *
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)    *
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF      *
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
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
