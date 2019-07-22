/********************************************************************
 * Copyright © 2019 Computational Molecular Biology Group,          *
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
 * « detailed description »
 *
 * @file MPIObservables.cpp
 * @brief « brief description »
 * @author chrisfroe
 * @date 03.06.19
 */

#include <readdy/kernel/mpi/observables/MPIObservables.h>
#include <readdy/kernel/mpi/MPIKernel.h>

namespace readdy::kernel::mpi::observables {

MPIVirial::MPIVirial(MPIKernel *kernel, Stride stride) : Virial(kernel, stride), kernel(kernel) {
}

void MPIVirial::evaluate() {
    // @todo MPI gather results, only save on master rank
    result = kernel->getMPIKernelStateModel().virial();
}

MPIPositions::MPIPositions(MPIKernel *kernel, unsigned int stride, const std::vector<std::string> &typesToCount)
        : Positions(kernel, stride, typesToCount), kernel(kernel) {}

void MPIPositions::evaluate() {
    // @todo MPI gather results, only save on master rank
    result.clear();
    auto &stateModel = kernel->getMPIKernelStateModel();
    const auto &pd = stateModel.getParticleData();
    auto it = pd->cbegin();
    if (typesToCount.empty()) {
        result = stateModel.getParticlePositions();
    } else {
        // only get positions of typesToCount
        while (it != pd->cend()) {
            if (!it->is_deactivated()) {
                if (std::find(typesToCount.begin(), typesToCount.end(), it->type) != typesToCount.end()) {
                    result.push_back(it->position());
                }
            }
            ++it;
        }
    }
}

MPIParticles::MPIParticles(MPIKernel *kernel, unsigned int stride) : Particles(kernel, stride), kernel(kernel) {}

void MPIParticles::evaluate() {
    // @todo MPI gather results, only save on master rank
    auto &resultTypes = std::get<0>(result);
    auto &resultIds = std::get<1>(result);
    auto &resultPositions = std::get<2>(result);
    resultTypes.clear();
    resultIds.clear();
    resultPositions.clear();
    const auto &particleData = kernel->getMPIKernelStateModel().getParticleData();
    auto it = particleData->cbegin();
    while (it != particleData->cend()) {
        if (!it->is_deactivated()) {
            resultTypes.push_back(it->type);
            resultIds.push_back(it->id);
            resultPositions.push_back(it->position());
        }
        ++it;
    }
}


MPIHistogramAlongAxis::MPIHistogramAlongAxis(MPIKernel *kernel, unsigned int stride,
                                             const std::vector<scalar> &binBorders,
                                             const std::vector<std::string> &typesToCount, unsigned int axis)
        : HistogramAlongAxis(kernel, stride, binBorders, typesToCount, axis), kernel(kernel) {}

void MPIHistogramAlongAxis::evaluate() {
    // @todo MPI gather results, only save on master rank
    std::fill(result.begin(), result.end(), 0);

    const auto &model = kernel->getMPIKernelStateModel();
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

MPINParticles::MPINParticles(MPIKernel *kernel, unsigned int stride, std::vector<std::string> typesToCount)
        : NParticles(kernel, stride, typesToCount), kernel(kernel) {}

void MPINParticles::evaluate() {
    // @todo MPI gather results, only save on master rank
    std::vector<unsigned long> resultVec = {};
    const auto &pd = kernel->getMPIKernelStateModel().getParticleData();

    if (typesToCount.empty()) {
        resultVec.push_back(pd->size() - pd->n_deactivated());
    } else {
        resultVec.resize(typesToCount.size());
        auto it = pd->cbegin();
        while (it != pd->cend()) {
            if (!it->is_deactivated()) {
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

MPIForces::MPIForces(MPIKernel *kernel, unsigned int stride, std::vector<std::string> typesToCount)
        : Forces(kernel, stride, typesToCount), kernel(kernel) {}

void MPIForces::evaluate() {
    // @todo MPI gather results, only save on master rank
    result.clear();
    const auto &pd = kernel->getMPIKernelStateModel().getParticleData();

    auto it = pd->cbegin();
    if (typesToCount.empty()) {
        // get all particles' forces
        for (; it != pd->cend(); ++it) {
            if (!it->is_deactivated()) {
                result.push_back(it->force);
            }
        }
    } else {
        // only get forces of typesToCount
        while (it != pd->cend()) {
            if (!it->is_deactivated()) {
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


MPIReactions::MPIReactions(MPIKernel *kernel, unsigned int stride) : Reactions(kernel, stride), kernel(kernel) {

}

void MPIReactions::evaluate() {
    // @todo MPI gather results, only save on master rank
    result = kernel->getMPIKernelStateModel().reactionRecords();
}

MPIReactionCounts::MPIReactionCounts(MPIKernel *kernel, unsigned int stride) : ReactionCounts(kernel, stride) {}

void MPIReactionCounts::evaluate() {
    // @todo MPI gather results, only save on master rank
    std::get<0>(result) = kernel->getMPIKernelStateModel().reactionCounts();

    // no topologies currently on MPI
    //std::get<1>(result) = kernel->getMPIKernelStateModel().spatialReactionCounts();
    //std::get<2>(result) = kernel->getMPIKernelStateModel().structuralReactionCounts();
}


}
