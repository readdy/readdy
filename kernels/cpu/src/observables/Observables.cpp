/**
 * << detailed description >>
 *
 * @file Observables.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 27.10.16
 */

#include <readdy/kernel/cpu/observables/Observables.h>
#include <readdy/kernel/cpu/CPUKernel.h>
#include <future>

namespace readdy {
namespace kernel {
namespace cpu {
namespace observables {

ParticlePosition::ParticlePosition(CPUKernel *const kernel, unsigned int stride,
                                   const std::vector<std::string> &typesToCount) :
        readdy::model::ParticlePositionObservable(kernel, stride, typesToCount), kernel(kernel) {}

void ParticlePosition::evaluate() {
    result.clear();
    const auto &pd = kernel->getKernelStateModel().getParticleData();
    const auto &entries = kernel->getKernelStateModel().getParticleData()->entries;
    if (typesToCount.empty()) {
        result = kernel->getKernelStateModel().getParticlePositions();
    } else {
        for (const auto &e : entries) {
            if (!e.is_deactivated() &&
                std::find(typesToCount.begin(), typesToCount.end(), e.type) != typesToCount.end()) {
                result.push_back(e.pos);
            }
        }
    }
}

HistogramAlongAxis::HistogramAlongAxis(CPUKernel *const kernel, unsigned int stride,
                                       const std::vector<double> &binBorders,
                                       const std::vector<std::string> &typesToCount, unsigned int axis)
        : readdy::model::HistogramAlongAxisObservable(kernel, stride, binBorders, typesToCount, axis),
          kernel(kernel) {
    size = result.size();
}

void HistogramAlongAxis::evaluate() {
    using Iter = readdy::kernel::cpu::model::ParticleData::entries_t::const_iterator;

    std::fill(result.begin(), result.end(), 0);

    const auto binBorders = this->binBorders;
    const auto typesToCount = this->typesToCount;
    const auto resultSize = result.size();
    const auto axis = this->axis;
    const auto data = kernel->getKernelStateModel().getParticleData();

    std::vector<std::future<result_t>> updates;
    auto worker = [binBorders, typesToCount, resultSize, data, axis](Iter from, Iter to, std::promise<result_t> update) {
        result_t resultUpdate;
        resultUpdate.resize(resultSize);

        const auto &entries = data->entries;

        for (auto it = from; it != to; ++it) {
            if (!it->is_deactivated() && typesToCount.find(it->type) != typesToCount.end()) {
                auto upperBound = std::upper_bound(binBorders.begin(), binBorders.end(), it->pos[axis]);
                if (upperBound != binBorders.end()) {
                    unsigned long binBordersIdx = upperBound - binBorders.begin();
                    if (binBordersIdx > 1) {
                        ++resultUpdate[binBordersIdx - 1];
                    }
                }
            }
        }

        update.set_value(std::move(resultUpdate));
    };

    {
        const std::size_t grainSize = size / kernel->getNThreads();

        std::vector<util::scoped_thread> threads;
        Iter workIter = data->entries.cbegin();
        for (unsigned int i = 0; i < kernel->getNThreads()-1; ++i) {
            std::promise<result_t> promise;
            updates.push_back(promise.get_future());
            threads.push_back(util::scoped_thread(std::thread(worker, workIter, workIter+grainSize, std::move(promise))));
            workIter+=grainSize;
        }
        std::promise<result_t> promise;
        updates.push_back(promise.get_future());
        threads.push_back(util::scoped_thread(std::thread(worker, workIter, data->entries.cend(), std::move(promise))));
    }

    for(auto& update : updates) {
        auto vec = std::move(update.get());
        std::transform(vec.begin(), vec.end(), result.begin(), result.end(), std::plus<result_t::value_type>());
    }
}


NParticles::NParticles(CPUKernel *const kernel, unsigned int stride, std::vector<std::string> typesToCount)
        : readdy::model::NParticlesObservable(kernel, stride, typesToCount),
          kernel(kernel) {}

void NParticles::evaluate() {
    std::vector<unsigned long> resultVec = {};
    if (typesToCount.empty()) {
        resultVec.push_back(kernel->getKernelStateModel().getParticleData()->size());
    } else {
        resultVec.resize(typesToCount.size());
        const auto &pd = kernel->getKernelStateModel().getParticleData();
        for (const auto &e : pd->entries) {
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

Forces::Forces(CPUKernel *const kernel, unsigned int stride, std::vector<std::string> typesToCount) :
        readdy::model::ForcesObservable(kernel, stride, typesToCount),
        kernel(kernel) {}

void Forces::evaluate() {
    result.clear();
    const auto &pd = kernel->getKernelStateModel().getParticleData();
    if (typesToCount.empty()) {
        result.reserve(pd->size());
    }
    for (const auto &e : pd->entries) {
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


}
}
}
}