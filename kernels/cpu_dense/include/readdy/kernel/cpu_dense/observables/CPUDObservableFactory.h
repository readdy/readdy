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


#include <readdy/model/observables/ObservableFactory.h>

namespace readdy {
namespace kernel {
namespace cpu_dense {
class CPUDKernel;
namespace observables {
class CPUDObservableFactory : public readdy::model::observables::ObservableFactory {

public:
    CPUDObservableFactory(CPUDKernel *const kernel);

    virtual readdy::model::observables::NParticles *
    createNParticlesObservable(unsigned int stride, std::vector<std::string> typesToCount = {}) const override;

    virtual readdy::model::observables::HistogramAlongAxis *
    createAxisHistogramObservable(unsigned int stride,
                                  std::vector<double> binBorders, std::vector<std::string> typesToCount,
                                  unsigned int axis) const override;

    virtual readdy::model::observables::Forces *
    createForcesObservable(unsigned int stride, std::vector<std::string> typesToCount = {}) const override;

    virtual readdy::model::observables::ParticlePosition *
    createParticlePositionObservable(unsigned int stride, std::vector<std::string> typesToCount = {}) const override;

private:
    CPUDKernel *const kernel;
};
}
}
}
}

#endif //READDY_DENSE_OBSERVABLEFACTORY_H
