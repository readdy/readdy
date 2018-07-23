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
