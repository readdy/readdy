/**
 * << detailed description >>
 *
 * @file ObservableFactory.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 23.11.16
 */


#include <readdy/kernel/cpu_dense/observables/CPUDObservableFactory.h>
#include <readdy/kernel/cpu_dense/CPUDKernel.h>
#include <readdy/kernel/cpu_dense/observables/CPUDObservables.h>

namespace readdy {
namespace kernel {
namespace cpu_dense {
namespace observables {
CPUDObservableFactory::CPUDObservableFactory(CPUDKernel *const kernel) : readdy::model::_internal::ObservableFactory(kernel),
                                                             kernel(kernel) {
}

readdy::model::NParticlesObservable *
CPUDObservableFactory::createNParticlesObservable(unsigned int stride, std::vector<std::string> typesToCount) const {
    return new CPUDNParticles(kernel, stride, typesToCount);
}

readdy::model::HistogramAlongAxisObservable *
CPUDObservableFactory::createAxisHistogramObservable(unsigned int stride, std::vector<double> binBorders,
                                                 std::vector<std::string> typesToCount,
                                                 unsigned int axis) const {
    return new CPUDHistogramAlongAxis(kernel, stride, binBorders, typesToCount, axis);
}

readdy::model::ForcesObservable *
CPUDObservableFactory::createForcesObservable(unsigned int stride, std::vector<std::string> typesToCount) const {
    return new CPUDForces(kernel, stride, typesToCount);
}

readdy::model::ParticlePositionObservable *
CPUDObservableFactory::createParticlePositionObservable(unsigned int stride, std::vector<std::string> typesToCount) const {
    return new CPUDParticlePosition(kernel, stride, typesToCount);
}

}
}
}
}