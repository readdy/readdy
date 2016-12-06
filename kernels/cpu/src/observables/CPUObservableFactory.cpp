/**
 * << detailed description >>
 *
 * @file ObservableFactory.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 21.07.16
 */

#include <readdy/kernel/cpu/observables/CPUObservableFactory.h>
#include <readdy/kernel/cpu/CPUKernel.h>
#include <readdy/kernel/cpu/observables/CPUObservables.h>
#include <readdy/kernel/singlecpu/observables/SCPUObservables.h>

namespace readdy {
namespace kernel {
namespace cpu {
namespace observables {
CPUObservableFactory::CPUObservableFactory(CPUKernel *const kernel) : readdy::model::observables::ObservableFactory(kernel),
                                                                kernel(kernel) {
}

readdy::model::observables::NParticles *
CPUObservableFactory::createNParticlesObservable(unsigned int stride, std::vector<std::string> typesToCount) const {
    return new CPUNParticles(kernel, stride, typesToCount);
}

readdy::model::observables::HistogramAlongAxis *
CPUObservableFactory::createAxisHistogramObservable(unsigned int stride, std::vector<double> binBorders,
                                                 std::vector<std::string> typesToCount,
                                                 unsigned int axis) const {
    return new CPUHistogramAlongAxis(kernel, stride, binBorders, typesToCount, axis);
}

readdy::model::observables::Forces *
CPUObservableFactory::createForcesObservable(unsigned int stride, std::vector<std::string> typesToCount) const {
    return new CPUForces(kernel, stride, typesToCount);
}

readdy::model::observables::ParticlePosition *
CPUObservableFactory::createParticlePositionObservable(unsigned int stride, std::vector<std::string> typesToCount) const {
    return new CPUParticlePosition(kernel, stride, typesToCount);
}

readdy::model::observables::RadialDistribution *
CPUObservableFactory::createRadialDistributionObservable(unsigned int stride, std::vector<double> binBorders, std::string typeCountFrom,
                                                      std::string typeCountTo, double particleToDensity) const {
    return new readdy::kernel::scpu::observables::RadialDistributionObservable<CPUKernel>(kernel, stride, binBorders, typeCountFrom, typeCountTo,
                                                                                               particleToDensity);
}

}
}
}
}
