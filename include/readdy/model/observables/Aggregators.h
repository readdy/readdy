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
 * An Aggregator is the same as a CombinerObservable but the results that it produces can accumulate and depend on time.
 * This file contains the definitions for aggregators, these are:
 *   - MeanSquaredDisplacement
 *
 * @file Aggregators.h
 * @brief Definition of several aggregators.
 * @author chrisfroe
 * @date 07.11.16
 */

#pragma once

#include <readdy/model/observables/Observable.h>
#include <readdy/model/observables/Observables.h>

NAMESPACE_BEGIN(readdy)
NAMESPACE_BEGIN(model)
NAMESPACE_BEGIN(observables)

class MeanSquaredDisplacement
        : public Combiner<std::pair<std::vector<time_step_type>, std::vector<scalar>>, Particles> {
public:
    MeanSquaredDisplacement(Kernel *kernel, stride_type stride, std::vector<std::string> typesToCount,
                            Particles *particlesObservable);

    MeanSquaredDisplacement(Kernel *kernel, stride_type stride, std::vector<ParticleTypeId> typesToCount,
                            Particles *particlesObservable);

    void evaluate() override = 0;

    std::string type() const override;

protected:
    std::vector<ParticleTypeId> typesToCount;
};

template<typename ParentObs>
class Trivial
        : public Combiner<std::pair<std::vector<time_step_type>, std::vector<typename ParentObs::result_type>>, ParentObs> {

    static_assert(
            std::is_base_of<readdy::model::observables::Observable<typename ParentObs::result_type>, ParentObs>::value,
            "ParentObs must extend readdy::model::observables::Observable");
public:

    using stride_type = typename ParentObs::stride_type;

    Trivial(Kernel *const kernel, stride_type stride, ParentObs *parentObs)
            : Combiner<std::pair<std::vector<time_step_type>, std::vector<typename ParentObs::result_type>>, ParentObs>(
            kernel, stride,
            parentObs) {}

    void evaluate() override {
        const auto &currentInput = std::get<0>(this->parentObservables)->getResult();
        auto &resultTimes = std::get<0>(this->result);
        auto &resultValues = std::get<1>(this->result);
        resultTimes.push_back(this->currentTimeStep());
        resultValues.push_back(currentInput);
    };

    std::string type() const override {
        return "Trivial";
    }
};

NAMESPACE_END(observables)
NAMESPACE_END(model)
NAMESPACE_END(readdy)
