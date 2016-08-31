/**
 * << detailed description >>
 *
 * @file ObservableFactory.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 21.07.16
 */

#include <readdy/kernel/cpu/observables/ObservableFactory.h>
#include <readdy/kernel/singlecpu/observables/SingleCPUObservables.h>
#include <readdy/kernel/cpu/CPUKernel.h>

namespace readdy {
    namespace kernel {
        namespace cpu {
            namespace observables {
                ObservableFactory::ObservableFactory(CPUKernel *const kernel) : readdy::model::_internal::ObservableFactory(kernel){
                }

                readdy::model::NParticlesObservable *
                ObservableFactory::createNParticlesObservable(readdy::model::Kernel *const kernel,
                                                              unsigned int stride, std::vector<std::string> typesToCount) const {
                    return new readdy::kernel::singlecpu::observables::NParticlesObservable<CPUKernel>(kernel, stride, typesToCount);
                }

                readdy::model::HistogramAlongAxisObservable *
                ObservableFactory::createAxisHistogramObservable(readdy::model::Kernel *const kernel,
                                                                 unsigned int stride, std::vector<double> binBorders,
                                                                 std::vector<std::string> typesToCount,
                                                                 unsigned int axis) const {
                    return new readdy::kernel::singlecpu::observables::SingleCPUHistogramAlongAxisObservable<readdy::kernel::cpu::CPUKernel>(kernel, stride, binBorders, typesToCount, axis);
                }

                readdy::model::ForcesObservable *
                ObservableFactory::createForcesObservable(readdy::model::Kernel *const kernel, unsigned int stride, std::vector<std::string> typesToCount) const {
                    return new readdy::kernel::singlecpu::observables::ForcesObservable<readdy::kernel::cpu::CPUKernel>(kernel, stride, typesToCount);
                }
            }
        }
    }
}