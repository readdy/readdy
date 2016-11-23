/**
 * << detailed description >>
 *
 * @file ObservableFactory.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 23.11.16
 */


#include <readdy/kernel/cpu_dense/observables/ObservableFactory.h>
#include <readdy/kernel/cpu_dense/Kernel.h>
#include <readdy/kernel/cpu_dense/observables/Observables.h>

namespace readdy {
namespace kernel {
namespace cpu_dense {
namespace observables {
ObservableFactory::ObservableFactory(Kernel *const kernel) : readdy::model::_internal::ObservableFactory(kernel),
                                                             kernel(kernel) {
}

readdy::model::NParticlesObservable *
ObservableFactory::createNParticlesObservable(unsigned int stride, std::vector<std::string> typesToCount) const {
    return new NParticles(kernel, stride, typesToCount);
}

readdy::model::HistogramAlongAxisObservable *
ObservableFactory::createAxisHistogramObservable(unsigned int stride, std::vector<double> binBorders,
                                                 std::vector<std::string> typesToCount,
                                                 unsigned int axis) const {
    return new HistogramAlongAxis(kernel, stride, binBorders, typesToCount, axis);
}

readdy::model::ForcesObservable *
ObservableFactory::createForcesObservable(unsigned int stride, std::vector<std::string> typesToCount) const {
    return new Forces(kernel, stride, typesToCount);
}

readdy::model::ParticlePositionObservable *
ObservableFactory::createParticlePositionObservable(unsigned int stride, std::vector<std::string> typesToCount) const {
    return new ParticlePosition(kernel, stride, typesToCount);
}

}
}
}
}