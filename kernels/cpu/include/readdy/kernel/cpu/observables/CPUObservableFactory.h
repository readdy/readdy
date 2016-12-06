/**
 * << detailed description >>
 *
 * @file ObservableFactory.h
 * @brief << brief description >>
 * @author clonker
 * @date 21.07.16
 */

#ifndef READDY_CPU_OBSERVABLEFACTORY_H
#define READDY_CPU_OBSERVABLEFACTORY_H


#include <readdy/model/observables/ObservableFactory.h>

namespace readdy {
namespace kernel {
namespace cpu {
class CPUKernel;
namespace observables {

class CPUObservableFactory : public readdy::model::observables::ObservableFactory {

public:
    CPUObservableFactory(CPUKernel *const kernel);

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

    virtual readdy::model::observables::RadialDistribution *
    createRadialDistributionObservable(unsigned int stride, std::vector<double> binBorders, std::string typeCountFrom,
                                        std::string typeCountTo, double particleToDensity) const override;

private:
    CPUKernel *const kernel;
};

}
}
}
}

#endif //READDY_CPU_OBSERVABLEFACTORY_H
