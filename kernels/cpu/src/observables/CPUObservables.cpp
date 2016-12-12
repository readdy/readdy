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

#include <readdy/common/thread/scoped_thread.h>

#include <readdy/kernel/cpu/observables/CPUObservables.h>
#include <readdy/kernel/cpu/CPUKernel.h>

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
    const auto &pd = kernel->getKernelStateModel().getParticleData();
    if (typesToCount.empty()) {
        result = kernel->getKernelStateModel().getParticlePositions();
    } else {
        for (const auto &e : *kernel->getKernelStateModel().getParticleData()) {
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
    const auto data = kernel->getKernelStateModel().getParticleData();

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

        std::vector<thd::scoped_thread> threads;
        Iter workIter = data->cbegin();
        for (unsigned int i = 0; i < kernel->getNThreads() - 1; ++i) {
            std::promise<result_t> promise;
            updates.push_back(promise.get_future());
            threads.push_back(thd::scoped_thread(std::thread(worker, workIter, workIter + grainSize, std::move(promise))));
            workIter += grainSize;
        }
        std::promise<result_t> promise;
        updates.push_back(promise.get_future());
        threads.push_back(thd::scoped_thread(std::thread(worker, workIter, data->cend(), std::move(promise))));
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
        resultVec.push_back(kernel->getKernelStateModel().getParticleData()->size());
    } else {
        resultVec.resize(typesToCount.size());
        const auto &pd = kernel->getKernelStateModel().getParticleData();
        for (const auto &e : *pd) {
            if (!e.is_deactivated()) {
                auto typeIt = std::find(typesToCount.begin(), typesToCount.end(), e.type);
                if (typeIt != typesToCount.end()) {
                    ++resultVec[typeIt - typesToCount.begin()];
                }
            }
        }
    }
    result = resultVec;
}

CPUForces::CPUForces(CPUKernel *const kernel, unsigned int stride, std::vector<std::string> typesToCount) :
        readdy::model::observables::Forces(kernel, stride, typesToCount),
        kernel(kernel) {}

void CPUForces::evaluate() {
    result.clear();
    const auto &pd = kernel->getKernelStateModel().getParticleData();
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
    const auto &particleData = kernel->getKernelStateModel().getParticleData();
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
}
}
}
}
