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
 * @file ObservableFactory.h
 * @brief Declaration of the observable factory for the CPU kernel
 * @author clonker
 * @date 21.07.16
 */

#pragma once

#include <readdy/model/observables/ObservableFactory.h>

namespace readdy::kernel::cpu {
class CPUKernel;
namespace observables {

class CPUObservableFactory : public readdy::model::observables::ObservableFactory {

public:
    explicit CPUObservableFactory(CPUKernel* kernel);

    [[nodiscard]] std::unique_ptr<model::observables::Energy>
    energy(Stride stride, ObsCallback <model::observables::Energy> callback) const override;

    [[nodiscard]] std::unique_ptr<model::observables::Virial>
    virial(Stride stride, ObsCallback <model::observables::Virial> callback) const override;

    [[nodiscard]] std::unique_ptr<model::observables::HistogramAlongAxis>
    histogramAlongAxis(Stride stride, std::vector<scalar> binBorders, std::vector<std::string> typesToCount,
                       unsigned int axis, ObsCallback <model::observables::HistogramAlongAxis> callback) const override;

    [[nodiscard]] std::unique_ptr<model::observables::NParticles>
    nParticles(Stride stride, std::vector<std::string> typesToCount,
               ObsCallback <model::observables::NParticles> callback) const override;

    [[nodiscard]] std::unique_ptr<model::observables::Forces>
    forces(Stride stride, std::vector<std::string> typesToCount,
           ObsCallback <model::observables::Forces> callback) const override;

    [[nodiscard]] std::unique_ptr<model::observables::Positions>
    positions(Stride stride, std::vector<std::string> typesToCount, ObsCallback<model::observables::Positions> callback) const override;

    [[nodiscard]] std::unique_ptr<model::observables::RadialDistribution>
    radialDistribution(Stride stride, std::vector<scalar> binBorders, std::vector<std::string> typeCountFrom,
                       std::vector<std::string> typeCountTo, scalar particleDensity,
                       ObsCallback<model::observables::RadialDistribution> callback) const override;

    [[nodiscard]] std::unique_ptr<model::observables::Particles>
    particles(Stride stride, ObsCallback<model::observables::Particles> callback) const override;

    [[nodiscard]] std::unique_ptr<model::observables::Reactions>
    reactions(Stride stride, ObsCallback<model::observables::Reactions> callback) const override;

    [[nodiscard]] std::unique_ptr<model::observables::ReactionCounts>
    reactionCounts(Stride stride, ObsCallback<model::observables::ReactionCounts> callback) const override;

private:
    CPUKernel *const kernel;
};

}
}
