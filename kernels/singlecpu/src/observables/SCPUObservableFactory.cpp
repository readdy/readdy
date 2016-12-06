/**
 * << detailed description >>
 *
 * @file SingleCPUObservableFactory.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 30.06.16
 */


#include <readdy/model/Kernel.h>
#include <readdy/kernel/singlecpu/observables/SCPUObservableFactory.h>
#include <readdy/kernel/singlecpu/observables/SCPUObservables.h>

namespace readdy {
namespace kernel {
namespace scpu {
namespace observables {
SCPUObservableFactory::SCPUObservableFactory(readdy::kernel::scpu::SCPUKernel *const kernel)
        : ObservableFactory(kernel), kernel(kernel) {
}

readdy::model::HistogramAlongAxisObservable *
SCPUObservableFactory::createAxisHistogramObservable(unsigned int stride, std::vector<double> binBorders,
                                                          std::vector<std::string> typesToCount,
                                                          unsigned int axis) const {
    return new SCPUHistogramAlongAxis(kernel, stride, binBorders, typesToCount, axis);
}

readdy::model::NParticlesObservable *SCPUObservableFactory::createNParticlesObservable(
        unsigned int stride, std::vector<std::string> typesToCount) const {
    return new SCPUNParticles(kernel, stride, typesToCount);
}

readdy::model::ForcesObservable *
SCPUObservableFactory::createForcesObservable(unsigned int stride, std::vector<std::string> typesToCount) const {
    return new SCPUForces(kernel, stride, typesToCount);
}

readdy::model::ParticlePositionObservable *
SCPUObservableFactory::createParticlePositionObservable(unsigned int stride, std::vector<std::string> typesToCount) const {
    return new SCPUParticlePosition(kernel, stride, typesToCount);
}

readdy::model::RadialDistributionObservable *
SCPUObservableFactory::createRadialDistributionObservable(unsigned int stride, std::vector<double> binBorders, std::string typeCountFrom,
                                                               std::string typeCountTo, double particleToDensity) const {
    return new RadialDistributionObservable<>(kernel, stride, binBorders, typeCountFrom, typeCountTo, particleToDensity);
}
}
}
}
}