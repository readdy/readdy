/**
 * << detailed description >>
 *
 * @file SingleCPUObservableFactory.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 30.06.16
 */


#include <readdy/model/Kernel.h>
#include <readdy/kernel/singlecpu/observables/SingleCPUObservableFactory.h>
#include <readdy/kernel/singlecpu/observables/SingleCPUObservables.h>

namespace readdy {
    namespace kernel {
        namespace singlecpu {
            namespace observables {
                SingleCPUObservableFactory::SingleCPUObservableFactory(readdy::model::Kernel *const kernel) : ObservableFactory(kernel) {
                }

                readdy::model::HistogramAlongAxisObservable *SingleCPUObservableFactory::createAxisHistogramObservable(readdy::model::Kernel *const kernel, unsigned int stride, std::vector<double> binBorders,
                                                                                                               std::vector<std::string> typesToCount, unsigned int axis) const {
                    return new SingleCPUHistogramAlongAxisObservable(kernel, stride, binBorders, typesToCount, axis);
                }


            }
        }
    }
}