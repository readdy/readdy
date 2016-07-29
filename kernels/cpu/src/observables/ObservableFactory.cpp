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
                                                              unsigned int stride) const {
                    return new readdy::kernel::singlecpu::observables::NParticlesObservable<CPUKernel>(kernel, stride);
                }
            }
        }
    }
}