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
 * @file ObservableFactory.h
 * @brief << brief description >>
 * @author clonker
 * @date 21.07.16
 */

#pragma once

#include <readdy/model/observables/ObservableFactory.h>

namespace readdy {
namespace kernel {
namespace cpu {
class CPUKernel;
namespace observables {

class CPUObservableFactory : public readdy::model::observables::ObservableFactory {

public:
    explicit CPUObservableFactory(CPUKernel* kernel);

    std::unique_ptr<model::observables::Virial> virial(stride_type stride) const override;

    std::unique_ptr<model::observables::HistogramAlongAxis>
    histogramAlongAxis(stride_type stride, std::vector<scalar> binBorders, std::vector<std::string> typesToCount,
                       unsigned int axis) const override;

    std::unique_ptr<model::observables::NParticles>
    nParticles(stride_type stride, std::vector<std::string> typesToCount) const override;

    std::unique_ptr<model::observables::Forces>
    forces(stride_type stride, std::vector<std::string> typesToCount) const override;

    std::unique_ptr<model::observables::Positions>
    positions(stride_type stride, std::vector<std::string> typesToCount) const override;

    std::unique_ptr<model::observables::RadialDistribution>
    radialDistribution(stride_type stride, std::vector<scalar> binBorders, std::vector<std::string> typeCountFrom,
                       std::vector<std::string> typeCountTo, scalar particleDensity) const override;

    std::unique_ptr<model::observables::Particles> particles(stride_type stride) const override;

    std::unique_ptr<model::observables::MeanSquaredDisplacement>
    msd(stride_type stride, std::vector<std::string> typesToCount,
        model::observables::Particles *particlesObservable) const override;

    std::unique_ptr<model::observables::Reactions> reactions(stride_type stride) const override;

    std::unique_ptr<model::observables::ReactionCounts> reactionCounts(stride_type stride) const override;

private:
    CPUKernel *const kernel;
};

}
}
}
}
