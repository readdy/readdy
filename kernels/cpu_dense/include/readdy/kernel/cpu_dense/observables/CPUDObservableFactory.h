/**
 * << detailed description >>
 *
 * @file ObservableFactory.h
 * @brief << brief description >>
 * @author clonker
 * @date 22.11.16
 */

#ifndef READDY_DENSE_OBSERVABLEFACTORY_H
#define READDY_DENSE_OBSERVABLEFACTORY_H


#include <readdy/model/_internal/ObservableFactory.h>

namespace readdy {
namespace kernel {
namespace cpu_dense {
class CPUDKernel;
namespace observables {
class CPUDObservableFactory : public readdy::model::_internal::ObservableFactory {

public:
    CPUDObservableFactory(CPUDKernel *const kernel);

    virtual readdy::model::NParticlesObservable *
    createNParticlesObservable(unsigned int stride, std::vector<std::string> typesToCount = {}) const override;

    virtual readdy::model::HistogramAlongAxisObservable *
    createAxisHistogramObservable(unsigned int stride,
                                  std::vector<double> binBorders, std::vector<std::string> typesToCount,
                                  unsigned int axis) const override;

    virtual readdy::model::ForcesObservable *
    createForcesObservable(unsigned int stride, std::vector<std::string> typesToCount = {}) const override;

    virtual readdy::model::ParticlePositionObservable *
    createParticlePositionObservable(unsigned int stride, std::vector<std::string> typesToCount = {}) const override;

private:
    CPUDKernel *const kernel;
};
}
}
}
}

#endif //READDY_DENSE_OBSERVABLEFACTORY_H
