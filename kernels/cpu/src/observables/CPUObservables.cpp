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
 * @file Observables.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 27.10.16
 */

#include <future>

#include <readdy/common/thread/scoped_async.h>

#include <readdy/kernel/cpu/observables/CPUObservables.h>
#include <readdy/kernel/cpu/CPUKernel.h>
#include <readdy/kernel/cpu/util/config.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace observables {

namespace thd = readdy::util::thread;

CPUPositions::CPUPositions(CPUKernel *const kernel, unsigned int stride,
                           const std::vector<std::string> &typesToCount) :
        readdy::model::observables::Positions(kernel, stride, typesToCount), kernel(kernel) {}

void CPUPositions::evaluate() {
    result.clear();
    auto &stateModel = kernel->getCPUKernelStateModel();
    const auto &pd = stateModel.getParticleData();
    if (typesToCount.empty()) {
        result = stateModel.getParticlePositions();
    } else {
        for (const auto &e : *stateModel.getParticleData()) {
            if (!e.is_deactivated() &&
                std::find(typesToCount.begin(), typesToCount.end(), e.type) != typesToCount.end()) {
                result.push_back(e.position());
            }
        }
    }
}

CPUHistogramAlongAxis::CPUHistogramAlongAxis(CPUKernel *const kernel, unsigned int stride,
                                             const std::vector<double> &binBorders,
                                             const std::vector<std::string> &typesToCount, unsigned int axis)
        : readdy::model::observables::HistogramAlongAxis(kernel, stride, binBorders, typesToCount, axis),
          kernel(kernel) {
    size = result.size();
}

void CPUHistogramAlongAxis::evaluate() {
    using Iter = readdy::kernel::cpu::model::CPUParticleData::entries_t::const_iterator;

    std::fill(result.begin(), result.end(), 0);

    const auto binBorders = this->binBorders;
    const auto typesToCount = this->typesToCount;
    const auto resultSize = result.size();
    const auto axis = this->axis;
    const auto data = kernel->getCPUKernelStateModel().getParticleData();

    std::vector<std::future<result_t>> updates;
    updates.reserve(kernel->getNThreads());
    auto worker = [binBorders, typesToCount, resultSize, data, axis](Iter from, Iter to, std::promise<result_t> update) {
        result_t resultUpdate;
        resultUpdate.resize(resultSize);

        for (Iter it = from; it != to; ++it) {
            if (!it->is_deactivated() && typesToCount.find(it->type) != typesToCount.end()) {
                auto upperBound = std::upper_bound(binBorders.begin(), binBorders.end(), it->position()[axis]);
                if (upperBound != binBorders.end()) {
                    unsigned long binBordersIdx = upperBound - binBorders.begin();
                    if (binBordersIdx >= 1 && binBordersIdx < resultSize) {
                        ++resultUpdate[binBordersIdx - 1];
                    }
                }
            }
        }

        update.set_value(std::move(resultUpdate));
    };

    {
        const std::size_t grainSize = data->size() / kernel->getNThreads();

        std::vector<threading_model> threads;
        Iter workIter = data->cbegin();
        for (unsigned int i = 0; i < kernel->getNThreads() - 1; ++i) {
            std::promise<result_t> promise;
            updates.push_back(promise.get_future());
            threads.emplace_back(worker, workIter, workIter + grainSize, std::move(promise));
            workIter += grainSize;
        }
        std::promise<result_t> promise;
        updates.push_back(promise.get_future());
        threads.emplace_back(worker, workIter, data->cend(), std::move(promise));
    }

    for (auto &update : updates) {
        auto vec = std::move(update.get());
        auto it1 = vec.begin();
        auto it2 = result.begin();
        for (; it1 != vec.end(); ++it1, ++it2) {
            *it2 += *it1;
        }
    }
}


CPUNParticles::CPUNParticles(CPUKernel *const kernel, unsigned int stride, std::vector<std::string> typesToCount)
        : readdy::model::observables::NParticles(kernel, stride, typesToCount),
          kernel(kernel) {}

void CPUNParticles::evaluate() {
    std::vector<unsigned long> resultVec = {};
    if (typesToCount.empty()) {
        resultVec.push_back(kernel->getCPUKernelStateModel().getParticleData()->size());
    } else {
        resultVec.resize(typesToCount.size());
        const auto &pd = kernel->getCPUKernelStateModel().getParticleData();
        for (const auto &e : *pd) {
            if (!e.is_deactivated()) {
                auto typeIt = std::find(typesToCount.begin(), typesToCount.end(), e.type);
                if (typeIt != typesToCount.end()) {
                    ++resultVec[typeIt - typesToCount.begin()];
                }
            }
        }
    }
    result = std::move(resultVec);
}

CPUForces::CPUForces(CPUKernel *const kernel, unsigned int stride, std::vector<std::string> typesToCount) :
        readdy::model::observables::Forces(kernel, stride, typesToCount),
        kernel(kernel) {}

void CPUForces::evaluate() {
    result.clear();
    const auto &pd = kernel->getCPUKernelStateModel().getParticleData();
    if (typesToCount.empty()) {
        result.reserve(pd->size());
    }
    for (const auto &e : *pd) {
        if (!e.is_deactivated()) {
            if (typesToCount.empty()) {
                result.push_back(e.force);
            } else {
                for (auto countedParticleType : typesToCount) {
                    if (e.type == countedParticleType) {
                        result.push_back(e.force);
                        break;
                    }
                }
            }
        }
    }
}


CPUParticles::CPUParticles(CPUKernel *const kernel, unsigned int stride)
        : readdy::model::observables::Particles(kernel, stride), kernel(kernel) {}

void CPUParticles::evaluate() {
    auto &resultTypes = std::get<0>(result);
    auto &resultIds = std::get<1>(result);
    auto &resultPositions = std::get<2>(result);
    resultTypes.clear();
    resultIds.clear();
    resultPositions.clear();
    const auto &particleData = kernel->getCPUKernelStateModel().getParticleData();
    resultTypes.reserve(particleData->size());
    resultIds.reserve(particleData->size());
    resultPositions.reserve(particleData->size());
    for (const auto &entry : *particleData) {
        if (!entry.is_deactivated()) {
            resultTypes.push_back(entry.type);
            resultIds.push_back(entry.id);
            resultPositions.push_back(entry.position());
        }
    }
}

CPUReactions::CPUReactions(CPUKernel *const kernel, unsigned int stride)
        : Reactions(kernel, stride), kernel(kernel) {}

void CPUReactions::evaluate() {
    const auto& model = kernel->getCPUKernelStateModel();
    const auto& records = model.reactionRecords();
    result.clear();
    result.reserve(records.size());
    result.insert(result.end(), records.begin(), records.end());
}

CPUReactionCounts::CPUReactionCounts(CPUKernel *const kernel, unsigned int stride)
        : ReactionCounts(kernel, stride), kernel(kernel) {}

void CPUReactionCounts::evaluate() {
    auto &order1 = std::get<0>(result);
    auto &order2 = std::get<1>(result);
    const auto& counts = kernel->getCPUKernelStateModel().reactionCounts();
    const auto& counts_order1 = std::get<0>(counts);
    const auto& counts_order2 = std::get<1>(counts);
    order1.assign(counts_order1.begin(), counts_order1.end());
    order2.assign(counts_order2.begin(), counts_order2.end());
}
}
}
}
}
