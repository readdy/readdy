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
 * @file Aggregators.cpp
 * @brief Definition of several aggregators.
 * @author chrisfroe
 * @date 15.11.16
 */

#include <readdy/model/observables/Aggregators.h>

#include <readdy/model/Kernel.h>


namespace readdy {
namespace model {
namespace observables {

MeanSquaredDisplacement::MeanSquaredDisplacement(Kernel *const kernel, stride_type stride,
                                                 std::vector<std::string> typesToCount,
                                                 Particles *particlesObservable)
        : MeanSquaredDisplacement(kernel, stride, readdy::model::_internal::util::transformTypes2(typesToCount,
                                                                                                  kernel->context()),
                                  particlesObservable) {}

MeanSquaredDisplacement::MeanSquaredDisplacement(Kernel *const kernel, stride_type stride,
                                                 std::vector<ParticleTypeId> typesToCount,
                                                 Particles *particlesObservable)
        : Combiner(kernel, stride, particlesObservable), typesToCount(std::move(typesToCount)) {}

std::string MeanSquaredDisplacement::type() const {
    return "Mean squared displacement";
}

}
}
}
