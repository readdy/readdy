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
    std::fill(result.begin(), result.end(), 0);

    const auto &model = kernel->getKernelStateModel();
    const auto data = model.getParticleData();
    const auto &entries = data->entries;

    for (const auto &e : entries) {
        if (!e.is_deactivated() && typesToCount.find(e.type) != typesToCount.end()) {
            const auto &vec = e.pos;
            auto upperBound = std::upper_bound(binBorders.begin(), binBorders.end(), vec[axis]);
            if (upperBound != binBorders.end()) {
                unsigned long binBordersIdx = upperBound - binBorders.begin();
                if (binBordersIdx > 1) {
                    ++result[binBordersIdx - 1];
                }
            }
        }
    }
}


NParticles::NParticles(CPUKernel *const kernel, unsigned int stride, std::vector<std::string> typesToCount)  : readdy::model::NParticlesObservable(kernel, stride, typesToCount),
                                                                                                               kernel(kernel) {}

void NParticles::evaluate()  {
    std::vector<unsigned long> resultVec = {};
    if (typesToCount.empty()) {
        resultVec.push_back(kernel->getKernelStateModel().getParticleData()->size());
    } else {
        resultVec.resize(typesToCount.size());
        const auto &pd = kernel->getKernelStateModel().getParticleData();
        for(const auto& e : pd->entries) {
            if(!e.is_deactivated()) {
                auto typeIt = std::find(typesToCount.begin(), typesToCount.end(), e.type);
                if(typeIt != typesToCount.end()) {
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
    if(typesToCount.empty()) {
        result.reserve(pd->size());
    }
    for(const auto& e : pd->entries) {
        if(!e.is_deactivated()) {
            if(typesToCount.empty()) {
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
