/**
 * << detailed description >>
 *
 * @file SingleCPUObservableFactory.h
 * @brief << brief description >>
 * @author clonker
 * @date 30.06.16
 */

#ifndef READDY_MAIN_SINGLECPUOBSERVABLEFACTORY_H
#define READDY_MAIN_SINGLECPUOBSERVABLEFACTORY_H

#include <readdy/model/_internal/ObservableFactory.h>
#include <readdy/kernel/singlecpu/SingleCPUKernel.h>

namespace readdy {
namespace kernel {
namespace singlecpu {
namespace observables {

class SingleCPUObservableFactory : public readdy::model::_internal::ObservableFactory {

public:
    SingleCPUObservableFactory(readdy::kernel::singlecpu::SingleCPUKernel *const kernel);

    virtual readdy::model::HistogramAlongAxisObservable *
    createAxisHistogramObservable(unsigned int stride, std::vector<double> binBorders,
                                  std::vector<std::string> typesToCount, unsigned int axis) const override;

    virtual readdy::model::NParticlesObservable *
    createNParticlesObservable(unsigned int stride, std::vector<std::string> typesToCount = {}) const override;

    virtual readdy::model::ForcesObservable *
    createForcesObservable(unsigned int stride, std::vector<std::string> typesToCount = {}) const override;

    virtual readdy::model::ParticlePositionObservable *
    createParticlePositionObservable(unsigned int stride, std::vector<std::string> typesToCount = {}) const override;

private:
    readdy::kernel::singlecpu::SingleCPUKernel *const kernel;
};

}
}
}
}
#endif //READDY_MAIN_SINGLECPUOBSERVABLEFACTORY_H
