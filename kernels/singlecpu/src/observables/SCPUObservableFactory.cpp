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
 * @file SCPUObservableFactory.cpp
 * @brief Implementation of observable factory of the single cpu kernel
 * @author clonker
 * @date 30.06.16
 */


#include <readdy/kernel/singlecpu/SCPUKernel.h>
#include <readdy/kernel/singlecpu/observables/SCPUObservableFactory.h>
#include <readdy/kernel/singlecpu/observables/SCPUObservables.h>

namespace readdy::kernel::scpu::observables {
SCPUObservableFactory::SCPUObservableFactory(readdy::kernel::scpu::SCPUKernel *const kernel)
        : ObservableFactory(kernel), kernel(kernel) {
}

std::unique_ptr<readdy::model::observables::HistogramAlongAxis>
SCPUObservableFactory::histogramAlongAxis(Stride stride, std::vector<scalar> binBorders,
                                          std::vector<std::string> typesToCount, unsigned int axis,
                                          ObsCallBack <readdy::model::observables::HistogramAlongAxis> callBack) const {
    auto obs = std::make_unique<SCPUHistogramAlongAxis>(kernel, stride, binBorders, typesToCount, axis);
    obs->setCallback(callBack);
    return std::move(obs);
}

std::unique_ptr<readdy::model::observables::NParticles>
SCPUObservableFactory::nParticles(Stride stride, std::vector<std::string> typesToCount,
                                  ObsCallBack <readdy::model::observables::NParticles> callback) const {
    auto obs = std::make_unique<SCPUNParticles>(kernel, stride, typesToCount);
    obs->setCallback(callback);
    return std::move(obs);
}

std::unique_ptr<readdy::model::observables::Forces>
SCPUObservableFactory::forces(Stride stride, std::vector<std::string> typesToCount,
                              ObsCallBack<readdy::model::observables::Forces> callback) const {
    auto obs = std::make_unique<SCPUForces>(kernel, stride, typesToCount);
    obs->setCallback(callback);
    return std::move(obs);
}

std::unique_ptr<readdy::model::observables::Positions>
SCPUObservableFactory::positions(Stride stride, std::vector<std::string> typesToCount,
                                 ObsCallBack<readdy::model::observables::Positions> callback) const {
    auto obs = std::make_unique<SCPUPositions>(kernel, stride, typesToCount);
    obs->setCallback(callback);
    return std::move(obs);
}

std::unique_ptr<readdy::model::observables::RadialDistribution>
SCPUObservableFactory::radialDistribution(Stride stride, std::vector<scalar> binBorders,
                                          std::vector<std::string> typeCountFrom,
                                          std::vector<std::string> typeCountTo, scalar particleDensity,
                                          ObsCallBack<readdy::model::observables::RadialDistribution> callback) const {
    auto obs = std::make_unique<readdy::model::observables::RadialDistribution>(
            kernel, stride, binBorders, typeCountFrom, typeCountTo, particleDensity);
    obs->setCallback(callback);
    return std::move(obs);
}

std::unique_ptr<readdy::model::observables::Particles>
SCPUObservableFactory::particles(Stride stride, ObsCallBack<readdy::model::observables::Particles> callback) const {
    auto obs = std::make_unique<SCPUParticles>(kernel, stride);
    obs->setCallback(callback);
    return std::move(obs);
}

std::unique_ptr<readdy::model::observables::Reactions>
SCPUObservableFactory::reactions(Stride stride, ObsCallBack<readdy::model::observables::Reactions> callback) const {
    auto obs = std::make_unique<SCPUReactions>(kernel, stride);
    obs->setCallback(callback);
    kernel->context().recordReactionsWithPositions() = true;
    return std::move(obs);
}

std::unique_ptr<readdy::model::observables::ReactionCounts>
SCPUObservableFactory::reactionCounts(Stride stride, ObsCallBack<readdy::model::observables::ReactionCounts> callback) const {
    auto obs = std::make_unique<SCPUReactionCounts>(kernel, stride);
    obs->setCallback(callback);
    kernel->context().recordReactionCounts() = true;
    return std::move(obs);
}

std::unique_ptr<readdy::model::observables::Virial>
SCPUObservableFactory::virial(Stride stride, ObsCallBack<readdy::model::observables::Virial> callBack) const {
    std::unique_ptr<readdy::model::observables::Virial> obs = std::make_unique<SCPUVirial>(kernel, stride);
    obs->setCallback(callBack);
    kernel->context().recordVirial() = true;
    return std::move(obs);
}

std::unique_ptr<readdy::model::observables::Energy> SCPUObservableFactory::energy(Stride stride, ObsCallBack<readdy::model::observables::Energy> callBack) const {
    auto obs = std::make_unique<readdy::model::observables::Energy>(kernel, stride);
    obs->setCallback(callBack);
    return std::move(obs);
}

}