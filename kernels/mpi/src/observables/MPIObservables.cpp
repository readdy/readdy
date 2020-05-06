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
 * Implementation of observables for the MPI kernel. In most cases during the evaluate(),
 * the workers collect the results which are then gathered on the master worker.
 *
 * @file MPIObservables.cpp
 * @brief Implementation of observables for the MPI kernel
 * @author chrisfroe
 * @date 03.06.19
 */

#include <utility>
#include <readdy/kernel/mpi/observables/MPIObservables.h>
#include <readdy/kernel/mpi/MPIKernel.h>

namespace readdy::kernel::mpi::observables {

MPIVirial::MPIVirial(MPIKernel *kernel, Stride stride) : Virial(kernel, stride), kernel(kernel) {}

void MPIVirial::evaluate() {
    // todo use MPI reduce with the standard '+' op, MPI_SUM
    //MPI_Reduce(send_data, recv_data, count, datatype, op, root, communicator);
    std::vector<readdy::model::observables::Virial::result_type> results(1);
    if (kernel->domain().isWorkerRank()) {
        results[0] = kernel->getMPIKernelStateModel().virial();
    }
    results = util::gatherObjects(results, 0, kernel->domain(), kernel->commUsedRanks());

    if (kernel->domain().isMasterRank()) {
        // add up virial tensors assuming that there was no double counting
        // which has to be ensured in calculateForces
        result = std::accumulate(results.begin(), results.end(), _internal::ReaDDyMatrix33<scalar>());
    }
}

void MPIVirial::append() {
    if (kernel->domain().isMasterRank()) {
        Virial::append();
    }
}

void MPIVirial::initializeDataSet(File &file, const std::string &dataSetName, Stride flushStride) {
    if (kernel->domain().isMasterRank()) {
        Virial::initializeDataSet(file, dataSetName, flushStride);
    }
}

MPIPositions::MPIPositions(MPIKernel *kernel, unsigned int stride, const std::vector<std::string> &typesToCount)
        : Positions(kernel, stride, typesToCount), kernel(kernel) {}

void MPIPositions::evaluate() {
    result.clear();
    if (kernel->domain().isWorkerRank()) {
        const auto &data = kernel->getMPIKernelStateModel().getParticleData();
        if (typesToCount.empty()) {
            for (const auto &p : *data) {
                if (!p.deactivated and p.responsible) {
                    result.push_back(p.pos);
                }
            }
        } else {
            // only get positions of typesToCount
            for (const auto &p : *data) {
                if (!p.deactivated and p.responsible) {
                    if (std::find(typesToCount.begin(), typesToCount.end(), p.type) != typesToCount.end()) {
                        result.push_back(p.pos);
                    }
                }
            }
        }
    }
    result = util::gatherObjects(result, 0, kernel->domain(), kernel->commUsedRanks());
}

void MPIPositions::append() {
    if (kernel->domain().isMasterRank()) {
        Positions::append();
    }
}

void MPIPositions::initializeDataSet(File &file, const std::string &dataSetName, Stride flushStride) {
    if (kernel->domain().isMasterRank()) {
        Positions::initializeDataSet(file, dataSetName, flushStride);
    }
}

MPIParticles::MPIParticles(MPIKernel *kernel, unsigned int stride) : Particles(kernel, stride), kernel(kernel) {}

void MPIParticles::evaluate() {
    auto &resultTypes = std::get<0>(result);
    auto &resultIds = std::get<1>(result);
    auto &resultPositions = std::get<2>(result);
    resultTypes.clear();
    resultIds.clear();
    resultPositions.clear();
    auto particles = kernel->getMPIKernelStateModel().gatherParticles();
    if (kernel->domain().isMasterRank()) {
        for (const auto &p : particles) {
            resultTypes.push_back(p.type());
            resultIds.push_back(p.id());
            resultPositions.push_back(p.pos());
        }
    }
}

void MPIParticles::append() {
    if (kernel->domain().isMasterRank()) {
        Particles::append();
    }
}

void MPIParticles::initializeDataSet(File &file, const std::string &dataSetName, Stride flushStride) {
    if (kernel->domain().isMasterRank()) {
        Particles::initializeDataSet(file, dataSetName, flushStride);
    }
}

MPIHistogramAlongAxis::MPIHistogramAlongAxis(MPIKernel *kernel, unsigned int stride,
                                             const std::vector<scalar> &binBorders,
                                             const std::vector<std::string> &typesToCount, unsigned int axis)
        : HistogramAlongAxis(kernel, stride, binBorders, typesToCount, axis), kernel(kernel) {}

void MPIHistogramAlongAxis::evaluate() {
    std::fill(result.begin(), result.end(), 0);
    if (kernel->domain().isWorkerRank()) {
        const auto data = kernel->getMPIKernelStateModel().getParticleData();
        for (const auto &p : *data) {
            if (!p.deactivated and p.responsible and typesToCount.find(p.type) != typesToCount.end()) {
                auto upperBound = std::upper_bound(binBorders.begin(), binBorders.end(), p.pos[axis]);
                if (upperBound != binBorders.end()) {
                    auto binBordersIdx = std::distance(binBorders.begin(), upperBound);
                    if (binBordersIdx > 1) {
                        ++result[binBordersIdx - 1];
                    }
                }
            }
        }
    }
    // todo variable float type
    HistogramAlongAxis::result_type tmp;
    tmp.resize(result.size());
    MPI_Reduce(result.data(), tmp.data(), static_cast<int>(result.size()), MPI_DOUBLE, MPI_SUM, 0, kernel->commUsedRanks());
    result = tmp;
}

void MPIHistogramAlongAxis::append() {
    if (kernel->domain().isMasterRank()) {
        HistogramAlongAxis::append();
    }
}

void MPIHistogramAlongAxis::initializeDataSet(File &file, const std::string &dataSetName, Stride flushStride) {
    if (kernel->domain().isMasterRank()) {
        HistogramAlongAxis::initializeDataSet(file, dataSetName, flushStride);
    }
}

MPINParticles::MPINParticles(MPIKernel *kernel, unsigned int stride, std::vector<std::string> typesToCount)
        : NParticles(kernel, stride, std::move(typesToCount)), kernel(kernel) {}

void MPINParticles::evaluate() {
    result.clear();
    if (kernel->domain().isWorkerRank()) {
        const auto &pd = kernel->getMPIKernelStateModel().getParticleData();
        if (typesToCount.empty()) {
            unsigned long n = std::count_if(pd->begin(), pd->end(),
                    [](const MPIEntry& entry)->bool{return (not entry.deactivated and entry.responsible);}
            );
            result.push_back(n);
        } else {
            result.resize(typesToCount.size());
            for (const auto &p : *pd) {
                if (!p.deactivated) {
                    unsigned int idx = 0;
                    for (const auto t : typesToCount) {
                        if (p.type == t) {
                            result[idx]++;
                            break;
                        }
                        ++idx;
                    }
                }
            }
        }
    } else if (kernel->domain().isMasterRank()) {
        if (typesToCount.empty()) {
            result.resize(1, 0);
        } else {
            result.resize(typesToCount.size(), 0);
        }
    } else {
        throw std::runtime_error("impossible");
    }

    NParticles::result_type tmp;
    tmp.resize(result.size());
    MPI_Reduce(result.data(), tmp.data(), static_cast<int>(result.size()), MPI_UNSIGNED_LONG, MPI_SUM, 0, kernel->commUsedRanks());
    result = tmp;
}

void MPINParticles::append() {
    if (kernel->domain().isMasterRank()) {
        NParticles::append();
    }
}

void MPINParticles::initializeDataSet(File &file, const std::string &dataSetName, Stride flushStride) {
    if (kernel->domain().isMasterRank()) {
        NParticles::initializeDataSet(file, dataSetName, flushStride);
    }
}

MPIForces::MPIForces(MPIKernel *kernel, unsigned int stride, std::vector<std::string> typesToCount)
        : Forces(kernel, stride, std::move(typesToCount)), kernel(kernel) {}

void MPIForces::evaluate() {
    result.clear();
    if (kernel->domain().isWorkerRank()) {
        const auto &pd = kernel->getMPIKernelStateModel().getParticleData();
        if (typesToCount.empty()) {
            // get all responsible particles' forces
            for (const auto &p : *pd) {
                if (!p.deactivated and p.responsible) {
                    result.push_back(p.force);
                }
            }
        } else {
            // only get forces of typesToCount
            for (const auto &p : *pd) {
                if (!p.deactivated and p.responsible) {
                    if (std::find(typesToCount.begin(), typesToCount.end(), p.type) != typesToCount.end()) {
                        result.push_back(p.force);
                    }
                }
            }
        }
    }
    result = util::gatherObjects(result, 0, kernel->domain(), kernel->commUsedRanks());
}

void MPIForces::append() {
    if (kernel->domain().isMasterRank()) {
        Forces::append();
    }
}

void MPIForces::initializeDataSet(File &file, const std::string &dataSetName, Stride flushStride) {
    if (kernel->domain().isMasterRank()) {
        Forces::initializeDataSet(file, dataSetName, flushStride);
    }
}

MPIReactions::MPIReactions(MPIKernel *kernel, unsigned int stride) : Reactions(kernel, stride), kernel(kernel) {}

void MPIReactions::evaluate() {
    result.clear();
    if (kernel->domain().isWorkerRank()) {
        result = kernel->getMPIKernelStateModel().reactionRecords();
    }
    result = util::gatherObjects(result, 0, kernel->domain(), kernel->commUsedRanks());
}

void MPIReactions::append() {
    if (kernel->domain().isMasterRank()) {
        Reactions::append();
    }
}

void MPIReactions::initializeDataSet(File &file, const std::string &dataSetName, Stride flushStride) {
    if (kernel->domain().isMasterRank()) {
        Reactions::initializeDataSet(file, dataSetName, flushStride);
    }
}

MPIReactionCounts::MPIReactionCounts(MPIKernel *kernel, unsigned int stride) : ReactionCounts(kernel, stride), kernel(kernel) {}

void MPIReactionCounts::evaluate() {
    // todo this could be nicer with MPI_Reduce: flatten into vec, mpi_reduce, copy into result map
    auto &counts = std::get<0>(result);

    // prepare a vector that mirrors the counts map
    std::vector<std::pair<readdy::ReactionId, unsigned long>> countsVec;
    if (kernel->domain().isWorkerRank()) {
        counts = kernel->getMPIKernelStateModel().reactionCounts();
        for (const auto &pair : counts) {
            countsVec.push_back(pair);
        }
    }

    // gather all vectors
    countsVec = util::gatherObjects(countsVec, 0, kernel->domain(), kernel->commUsedRanks());

    // reduce the countsVec on master rank, i.e. aggregate pairs with the same id into the usual map format
    if (kernel->domain().isMasterRank()) {
        counts.clear();
        for (const auto &[id, n] : countsVec) {
            if (counts.find(id) != counts.end()) {
                counts[id] += n;
            } else {
                counts[id] = n;
            }
        }
    }

    // no topologies currently on MPI
    //std::get<1>(result) = kernel->getMPIKernelStateModel().spatialReactionCounts();
    //std::get<2>(result) = kernel->getMPIKernelStateModel().structuralReactionCounts();
}

void MPIReactionCounts::append() {
    if (kernel->domain().isMasterRank()) {
        ReactionCounts::append();
    }
}

void MPIReactionCounts::initializeDataSet(File &file, const std::string &dataSetName, Stride flushStride) {
    if (kernel->domain().isMasterRank()) {
        ReactionCounts::initializeDataSet(file, dataSetName, flushStride);
    }
}

MPIEnergy::MPIEnergy(MPIKernel *kernel, Stride stride) : Energy(kernel, stride), kernel(kernel) {}

void MPIEnergy::evaluate() {
    result = 0.;
    if (kernel->domain().isWorkerRank()) {
        result = kernel->stateModel().energy();
    }
    readdy::scalar tmp;
    MPI_Reduce(&result, &tmp, static_cast<int>(1), MPI_DOUBLE, MPI_SUM, 0, kernel->commUsedRanks());
    result = tmp;
}

void MPIEnergy::append() {
    if (kernel->domain().isMasterRank()) {
        Energy::append();
    }
}

void MPIEnergy::initializeDataSet(File &file, const std::string &dataSetName, Stride flushStride) {
    if (kernel->domain().isMasterRank()) {
        Energy::initializeDataSet(file, dataSetName, flushStride);
    }
}

}
