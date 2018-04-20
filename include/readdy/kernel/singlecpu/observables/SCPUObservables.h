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
 * @file SingleCPUObservables.h
 * @brief << brief description >>
 * @author clonker
 * @date 30.06.16
 */

#pragma once

#include <readdy/model/observables/Observables.h>
#include <readdy/kernel/singlecpu/SCPUKernel.h>


namespace readdy {
namespace kernel {
namespace scpu {
namespace observables {

class SCPUPositions : public readdy::model::observables::Positions {
public:
    SCPUPositions(SCPUKernel *const kernel, unsigned int stride, const std::vector<std::string> &typesToCount = {}) :
            readdy::model::observables::Positions(kernel, stride, typesToCount), kernel(kernel) {}

    void evaluate() override {
        result.clear();
        auto& stateModel = kernel->getSCPUKernelStateModel();
        const auto &pd = stateModel.getParticleData();
        auto it = pd->cbegin();
        if (typesToCount.empty()) {
            result = stateModel.getParticlePositions();
        } else {
            // only get positions of typesToCount
            while (it != pd->cend()) {
                if(!it->is_deactivated()) {
                    if (std::find(typesToCount.begin(), typesToCount.end(), it->type) != typesToCount.end()) {
                        result.push_back(it->position());
                    }
                }
                ++it;
            }
        }
    }

protected:
    SCPUKernel *const kernel;
};

class SCPUParticles : public readdy::model::observables::Particles {
public:
    SCPUParticles(SCPUKernel *const kernel, unsigned int stride)
            : readdy::model::observables::Particles(kernel, stride), kernel(kernel) {};

    void evaluate() override {
        auto &resultTypes = std::get<0>(result);
        auto &resultIds = std::get<1>(result);
        auto &resultPositions = std::get<2>(result);
        resultTypes.clear();
        resultIds.clear();
        resultPositions.clear();
        const auto &particleData = kernel->getSCPUKernelStateModel().getParticleData();
        auto it = particleData->cbegin();
        while(it != particleData->cend()) {
            if(!it->is_deactivated()) {
                resultTypes.push_back(it->type);
                resultIds.push_back(it->id);
                resultPositions.push_back(it->position());
            }
            ++it;
        }
    };

protected:
    SCPUKernel *const kernel;
};

class SCPUHistogramAlongAxis : public readdy::model::observables::HistogramAlongAxis {

public:
    SCPUHistogramAlongAxis(SCPUKernel *const kernel, unsigned int stride,
                           const std::vector<scalar> &binBorders,
                           const std::vector<std::string> &typesToCount,
                           unsigned int axis)
            : readdy::model::observables::HistogramAlongAxis(kernel, stride, binBorders, typesToCount, axis),
              kernel(kernel) {
        size = result.size();
    }

    void evaluate() override {
        std::fill(result.begin(), result.end(), 0);

        const auto &model = kernel->getSCPUKernelStateModel();
        const auto data = model.getParticleData();

        auto it = data->cbegin();

        while (it != data->cend()) {
            if (!it->is_deactivated() and typesToCount.find(it->type) != typesToCount.end()) {
                const auto &vec = it->position();
                auto upperBound = std::upper_bound(binBorders.begin(), binBorders.end(), vec[axis]);
                if (upperBound != binBorders.end()) {
                    unsigned long binBordersIdx = static_cast<unsigned long>(upperBound - binBorders.begin());
                    if (binBordersIdx > 1) {
                        ++result[binBordersIdx - 1];
                    }
                }
            }
            ++it;
        }
    }

protected:
    SCPUKernel *const kernel;
    size_t size;
};

class SCPUNParticles : public readdy::model::observables::NParticles {
public:
    SCPUNParticles(SCPUKernel *const kernel, unsigned int stride, const std::vector<std::string>& typesToCount = {}) :
            readdy::model::observables::NParticles(kernel, stride, typesToCount),
            singleCPUKernel(kernel) {}

    void evaluate() override {
        std::vector<unsigned long> resultVec = {};
        const auto &pd = singleCPUKernel->getSCPUKernelStateModel().getParticleData();

        if (typesToCount.empty()) {
            resultVec.push_back(pd->size() - pd->n_deactivated());
        } else {
            resultVec.resize(typesToCount.size());
            auto it = pd->cbegin();
            while(it != pd->cend()) {
                if(!it->is_deactivated()) {
                    unsigned int idx = 0;
                    for (const auto t : typesToCount) {
                        if (it->type == t) {
                            resultVec[idx]++;
                            break;
                        }
                        ++idx;
                    }
                }
                ++it;
            }
        }
        result = resultVec;
    }

protected:
    SCPUKernel *const singleCPUKernel;
};

class SCPUForces : public readdy::model::observables::Forces {
public:
    SCPUForces(SCPUKernel *const kernel, unsigned int stride, const std::vector<std::string> &typesToCount = {}) :
            readdy::model::observables::Forces(kernel, stride, typesToCount),
            kernel(kernel) {}

    ~SCPUForces() override = default;
    SCPUForces(const SCPUForces&) = default;
    SCPUForces& operator=(const SCPUForces&) = default;
    SCPUForces(SCPUForces&&) = default;
    SCPUForces& operator=(SCPUForces&&) = default;

    void evaluate() override {
        result.clear();
        const auto &pd = kernel->getSCPUKernelStateModel().getParticleData();

        auto it = pd->cbegin();
        if (typesToCount.empty()) {
            // get all particles' forces
            for(; it != pd->cend(); ++it) {
                if(!it->is_deactivated()) {
                    result.push_back(it->force);
                }
            }
        } else {
            // only get forces of typesToCount
            while (it != pd->cend()) {
                if(!it->is_deactivated()) {
                    for (auto countedParticleType : typesToCount) {
                        if (it->type == countedParticleType) {
                            result.push_back(it->force);
                            break;
                        }
                    }
                }
                ++it;
            }
        }
    }


protected:
    SCPUKernel *const kernel;
};

class SCPUVirial : public readdy::model::observables::Virial {
public:
    SCPUVirial(SCPUKernel *kernel, stride_type stride) : Virial(kernel, stride), _kernel(kernel) {};

    void evaluate() override {
        result = _kernel->getSCPUKernelStateModel().virial();
    }

private:
    SCPUKernel* const _kernel;
};

class SCPUReactions : public readdy::model::observables::Reactions {
public:
    SCPUReactions(SCPUKernel *const kernel, unsigned int stride)
            : Reactions(kernel, stride), kernel(kernel) {}

    ~SCPUReactions() override = default;

    void evaluate() override {
        const auto &records = kernel->getSCPUKernelStateModel().reactionRecords();
        result = records;
    }

private:
    SCPUKernel* const kernel;
};

class SCPUReactionCounts : public readdy::model::observables::ReactionCounts {
public:
    SCPUReactionCounts(SCPUKernel *const kernel, unsigned int stride) : ReactionCounts(kernel, stride), kernel(kernel) {}

    ~SCPUReactionCounts() override = default;

    void evaluate() override {
        result = kernel->getSCPUKernelStateModel().reactionCounts();
    };

private:
    SCPUKernel* const kernel;
};

}
}
}
}
