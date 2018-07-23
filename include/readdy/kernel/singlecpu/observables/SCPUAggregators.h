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
 * @file Aggregators.h
 * @brief SingleCpu implementations of the aggregators defined in the core library.
 * @author chrisfroe
 * @date 09.11.16
 */

#pragma once
#include <readdy/model/observables/Aggregators.h>
#include <readdy/kernel/singlecpu/SCPUKernel.h>

namespace readdy {
namespace kernel {
namespace scpu {
namespace observables {

// todo check what happens with ids if particle reacts
/**
 * MSD calculation, that remembers particle ids and initial positions on first call. In subsequent calls it only searches
 * for those particle ids. I.e. if all initial particles vanished, there is nothing to compute anymore.
 */

template<typename KERNEL=readdy::kernel::scpu::SCPUKernel>
class SCPUMeanSquaredDisplacement : public readdy::model::observables::MeanSquaredDisplacement {
public:
    SCPUMeanSquaredDisplacement(KERNEL *const kernel, unsigned int stride, std::vector<std::string> typesToCount,
                            readdy::model::observables::Particles *particlesObservable)
            : readdy::model::observables::MeanSquaredDisplacement(kernel, stride, typesToCount, particlesObservable), kernel(kernel) {};

    void evaluate() override {
        const auto &currentInput = std::get<0>(parentObservables)->getResult();
        const auto &types = std::get<0>(currentInput);
        const auto &ids = std::get<1>(currentInput);
        const auto &positions = std::get<2>(currentInput);

        auto &resultTime = std::get<0>(result);
        resultTime.push_back(currentTimeStep());
        auto &resultMsd = std::get<1>(result);
        if (resultMsd.size() == 0) { // first evaluation
            for (auto i = 0; i < positions.size(); ++i) {
                if (typesToCount.size() == 0 || std::find(typesToCount.begin(), typesToCount.end(), types[i]) != typesToCount.end()) {
                    initialPositions.emplace(std::make_pair(ids[i], positions[i]));
                }
            }
            resultMsd.push_back(0);
            numberOfParticles.push_back(initialPositions.size());
            return;
        }
        scalar msd = 0;
        unsigned long currentNumberParticles = 0;
        auto posIt = positions.begin();
        auto idsIt = ids.begin();
        for(; idsIt != ids.end(); ++idsIt, ++posIt) {
            auto it = initialPositions.find(*idsIt);
            if (it != initialPositions.end()) {
                const auto displacement = *posIt - it->second;
                msd += displacement * displacement;
                ++currentNumberParticles;
            }
        }
        if (currentNumberParticles > 0) {
            msd /= static_cast<scalar>(currentNumberParticles);
        };
        resultMsd.push_back(msd);
        numberOfParticles.push_back(currentNumberParticles);
    }

protected:
    std::unordered_map<readdy::model::Particle::id_type, Vec3> initialPositions;
    std::vector<unsigned long> numberOfParticles;
    KERNEL *const kernel;
};

}
}
}
}
