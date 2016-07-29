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
#include <readdy/kernel/singlecpu/SingleCPUKernel.h>


namespace readdy {
    namespace kernel {
        namespace singlecpu {
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

                template<typename kernel_t=readdy::kernel::singlecpu::SingleCPUKernel>
                class NParticlesObservable : public readdy::model::NParticlesObservable {
                public:
                    NParticlesObservable(readdy::model::Kernel *const kernel, unsigned int stride) :
                            readdy::model::NParticlesObservable(kernel, stride),
                            singleCPUKernel(dynamic_cast<kernel_t *>(kernel))
                    { }

                    virtual void evaluate() override {
                        result = singleCPUKernel->getKernelStateModel().getParticleData()->size();
                    }

                protected:
                     kernel_t*const singleCPUKernel;
                };

            }
        }
    }
}
#endif //READDY_MAIN_SINGLECPUOBSERVABLES_H
