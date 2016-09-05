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


#include <readdy/model/_internal/ObservableFactory.h>

namespace readdy {
    namespace kernel {
        namespace cpu {
            class CPUKernel;
            namespace observables {
                class ObservableFactory : public readdy::model::_internal::ObservableFactory {

                public:
                    ObservableFactory(CPUKernel *const kernel);
                    virtual readdy::model::NParticlesObservable *
                    createNParticlesObservable(unsigned int stride, std::vector<std::string> typesToCount = {}) const override;

                    virtual readdy::model::HistogramAlongAxisObservable *
                    createAxisHistogramObservable(unsigned int stride,
                                                  std::vector<double> binBorders, std::vector<std::string> typesToCount,
                                                  unsigned int axis) const override;

                    virtual readdy::model::ForcesObservable *
                    createForcesObservable(unsigned int stride, std::vector<std::string> typesToCount = {}) const override;

                private:
                    CPUKernel *const kernel;
                };
            }
        }
    }
}

#endif //READDY_CPU_OBSERVABLEFACTORY_H
