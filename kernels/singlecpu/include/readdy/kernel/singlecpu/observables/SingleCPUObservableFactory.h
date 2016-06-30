/**
 * << detailed description >>
 *
 * @file SingleCPUObservableFactory.h
 * @brief << brief description >>
 * @author clonker
 * @date 30.06.16
 */

#ifndef READDY_MAIN_SINGLECPUOBSERVABLEFACTORY_H
#define READDY_MAIN_SINGLECPUOBSERVABLEFACTORY_H

#include <readdy/model/_internal/ObservableFactory.h>

namespace readdy {
    namespace kernel {
        namespace singlecpu {
            namespace observables {
                class SingleCPUObservableFactory : public readdy::model::_internal::ObservableFactory{

                public:
                    SingleCPUObservableFactory(readdy::model::Kernel *const kernel);

                    virtual readdy::model::HistogramAlongAxisObservable *createAxisHistogramObservable(readdy::model::Kernel *const kernel, unsigned int stride, std::vector<double> binBorders,
                                                                                               std::vector<std::string> typesToCount, unsigned int axis) const override;

                };
            }
        }
    }
}
#endif //READDY_MAIN_SINGLECPUOBSERVABLEFACTORY_H
