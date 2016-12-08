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

#ifndef READDY_MAIN_SCPUAGGREGATORS_H
#define READDY_MAIN_SCPUAGGREGATORS_H

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

class SCPUMeanSquaredDisplacement : public readdy::model::observables::MeanSquaredDisplacement {
public:
    SCPUMeanSquaredDisplacement(SCPUKernel *const kernel, unsigned int stride, std::vector<std::string> typesToCount,
                            readdy::model::observables::Particles *particlesObservable)
            : readdy::model::observables::MeanSquaredDisplacement(kernel, stride, typesToCount, particlesObservable), kernel(kernel) {};

    virtual void evaluate() override {
        const auto &currentInput = std::get<0>(parentObservables)->getResult();
        const auto &types = std::get<0>(currentInput);
        const auto &ids = std::get<1>(currentInput);
        const auto &positions = std::get<2>(currentInput);

        auto &resultTime = std::get<0>(result);
        resultTime.push_back(getCurrentTimeStep());
        auto &resultMsd = std::get<1>(result);
        if (resultMsd.size() == 0) { // first evaluation
            for (auto i = 0; i < positions.size(); ++i) {
                if (typesToCount.size() == 0 || std::find(typesToCount.begin(), typesToCount.end(), types[i]) != typesToCount.end()) {
                    initialPositions.emplace(std::make_pair(ids[i], positions[i]));
                }
            }
            resultMsd.push_back(0);
            return;
        }
        double msd = 0;
        unsigned int numberParticles = 0;
        auto initPos = initialPositions.begin();
        while (initPos != initialPositions.end()) {
            auto it = std::find(ids.begin(), ids.end(), (*initPos).first);
            if (it != ids.end()) {
                // idx is the position in the current(!) result of the particle with a certain id
                auto idx = it - ids.begin();
                const auto displacement = positions[idx] - (*initPos).second;
                msd += displacement * displacement;
                ++numberParticles;
                ++initPos;
            } else {
                // particle does not exist anymore -> remove it from map and DO NOT increment the iterator
                initPos = initialPositions.erase(initPos);
            }
        }
        if (numberParticles == 0) numberParticles = 1;
        msd /= static_cast<double>(numberParticles);
        resultMsd.push_back(msd);
    }

protected:
    std::unordered_map<readdy::model::Particle::id_type, readdy::model::Vec3> initialPositions;
    SCPUKernel *const kernel;
};

}
}
}
}

#endif //READDY_MAIN_SCPUAGGREGATORS_H
