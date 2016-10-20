/**
 * << detailed description >>
 *
 * @file SingleCPUObservables.h
 * @brief << brief description >>
 * @author clonker
 * @date 30.06.16
 */


#ifndef READDY_MAIN_SINGLECPUOBSERVABLES_H
#define READDY_MAIN_SINGLECPUOBSERVABLES_H

#include <readdy/model/Observables.h>
#include <readdy/kernel/singlecpu/SingleCPUKernel.h>


namespace readdy {
namespace kernel {
namespace singlecpu {
namespace observables {

template<typename kernel_t=readdy::kernel::singlecpu::SingleCPUKernel>
class ParticlePositionObservable : public readdy::model::ParticlePositionObservable {
public:
    ParticlePositionObservable(kernel_t *const kernel, unsigned int stride, const std::vector<std::string> &typesToCount = {}) :
            readdy::model::ParticlePositionObservable(kernel, stride, typesToCount), kernel(kernel) {}

    virtual void evaluate() override {
        result.clear();
        const auto &pd = kernel->getKernelStateModel().getParticleData();
        auto positionsIt = pd->cbegin_positions();
        if (typesToCount.empty()) {
            // get all particles' positions
            result.reserve(pd->size());
            std::copy(positionsIt, pd->cend_positions(), std::back_inserter(result));
        } else {
            // only get positions of typesToCount
            auto typesIt = pd->cbegin_types();
            while (positionsIt != pd->cend_positions()) {
                if (std::find(typesToCount.begin(), typesToCount.end(), *typesIt) != typesToCount.end()) {
                    result.push_back(*positionsIt);
                }
                ++positionsIt;
                ++typesIt;
            }
        }
    }

protected:
    kernel_t *const kernel;
};

template<typename kernel_t=readdy::kernel::singlecpu::SingleCPUKernel>
class HistogramAlongAxisObservable : public readdy::model::HistogramAlongAxisObservable {

public:
    HistogramAlongAxisObservable(kernel_t *const kernel, unsigned int stride,
                                 const std::vector<double> &binBorders,
                                 const std::vector<std::string> &typesToCount,
                                 unsigned int axis)
            : readdy::model::HistogramAlongAxisObservable(kernel, stride, binBorders, typesToCount, axis),
              kernel(kernel) {
        size = result.size();
    }

    virtual void evaluate() override {
        std::fill(result.begin(), result.end(), 0);

        const auto &model = kernel->getKernelStateModel();
        const auto data = model.getParticleData();

        auto it_pos = data->begin_positions();
        auto it_types = data->begin_types();

        while (it_pos != data->end_positions()) {
            if (typesToCount.find(*it_types) != typesToCount.end()) {
                const auto &vec = *it_pos;
                auto upperBound = std::upper_bound(binBorders.begin(), binBorders.end(), vec[axis]);
                if (upperBound != binBorders.end()) {
                    unsigned long binBordersIdx = upperBound - binBorders.begin();
                    if (binBordersIdx > 1) {
                        ++result[binBordersIdx - 1];
                    }
                }
            }
            ++it_pos;
            ++it_types;
        }
    }

protected:
    kernel_t *const kernel;
    size_t size;
};

template<typename kernel_t=readdy::kernel::singlecpu::SingleCPUKernel>
class NParticlesObservable : public readdy::model::NParticlesObservable {
public:
    NParticlesObservable(kernel_t *const kernel, unsigned int stride, std::vector<std::string> typesToCount = {}) :
            readdy::model::NParticlesObservable(kernel, stride, typesToCount),
            singleCPUKernel(kernel) {}

    virtual void evaluate() override {
        std::vector<unsigned long> resultVec = {};
        if (typesToCount.empty()) {
            resultVec.push_back(singleCPUKernel->getKernelStateModel().getParticleData()->size());
        } else {
            resultVec.resize(typesToCount.size());
            const auto &pd = singleCPUKernel->getKernelStateModel().getParticleData();
            auto typesIt = pd->cbegin_types();
            while (typesIt != pd->cend_types()) {
                unsigned int idx = 0;
                for (const auto t : typesToCount) {
                    if (*typesIt == t) {
                        resultVec[idx]++;
                        break;
                    }
                    ++idx;
                }
                ++typesIt;
            }
        }
        result = resultVec;
    }

protected:
    kernel_t *const singleCPUKernel;
};

template<typename kernel_t=readdy::kernel::singlecpu::SingleCPUKernel>
class ForcesObservable : public readdy::model::ForcesObservable {
public:
    ForcesObservable(kernel_t *const kernel, unsigned int stride, std::vector<std::string> typesToCount = {}) :
            readdy::model::ForcesObservable(kernel, stride, typesToCount),
            kernel(kernel) {}

    virtual ~ForcesObservable() {}

    virtual void evaluate() override {
        result.clear();
        const auto &pd = kernel->getKernelStateModel().getParticleData();
        auto forcesIt = pd->cbegin_forces();

        if (typesToCount.empty()) {
            // get all particles' forces
            result.reserve(pd->size());
            std::copy(forcesIt, pd->cend_forces(), std::back_inserter(result));
        } else {
            // only get forces of typesToCount
            auto typesIt = pd->cbegin_types();
            while (forcesIt != pd->cend_forces()) {
                for (auto countedParticleType : typesToCount) {
                    if (*typesIt == countedParticleType) {
                        result.push_back(*forcesIt);
                        break;
                    }
                }
                ++forcesIt;
                ++typesIt;
            }
        }
    }


protected:
    kernel_t *const kernel;
};

}
}
}
}
#endif //READDY_MAIN_SINGLECPUOBSERVABLES_H
