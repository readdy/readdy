/**
 * << detailed description >>
 *
 * @file SingleCPUObservables.h
 * @brief << brief description >>
 * @author clonker
 * @date 30.06.16
 */


#ifndef READDY_MAIN_SINGLECPUOBSERVABLES_H
#define READDY_MAIN_SINGLECPUOBSERVABLES_H

#include <readdy/model/Observables.h>


namespace readdy {
    namespace kernel {
        namespace singlecpu {
            class SingleCPUKernel;
            namespace observables {
                class SingleCPUHistogramAlongAxisObservable : public readdy::model::HistogramAlongAxisObservable {

                public:
                    SingleCPUHistogramAlongAxisObservable(readdy::model::Kernel *const kernel, unsigned int stride, const std::vector<double> &binBorders, const std::vector<std::string> &typesToCount,
                                                          unsigned int axis);

                    virtual void evaluate() override;

                protected:
                    readdy::kernel::singlecpu::SingleCPUKernel *const singleCPUKernel;
                    size_t size;
                };

                class SingleCPUNparticlesObservable : public readdy::model::NParticlesObservable {

                public:
                    SingleCPUNparticlesObservable(readdy::model::Kernel *const kernel, unsigned int stride);

                    virtual void evaluate() override;

                protected:
                    readdy::kernel::singlecpu::SingleCPUKernel *const singleCPUKernel;

                };

            }
        }
    }
}
#endif //READDY_MAIN_SINGLECPUOBSERVABLES_H
