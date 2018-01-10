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
