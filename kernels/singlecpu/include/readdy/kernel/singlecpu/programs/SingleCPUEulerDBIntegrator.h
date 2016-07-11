/**
 * << detailed description >>
 *
 * @file SingleCPUEulerDBIntegrator.h
 * @brief << brief description >>
 * @author clonker
 * @date 19.04.16
 */

#ifndef READDY_MAIN_SingleCPUEulerDBIntegrator_H
#define READDY_MAIN_SingleCPUEulerDBIntegrator_H

#include <readdy/model/programs/Programs.h>
#include <readdy/kernel/singlecpu/SingleCPUKernel.h>

namespace readdy {
    namespace kernel {
        namespace singlecpu {

            namespace programs {
                class SingleCPUEulerDBIntegrator : public readdy::model::programs::EulerBDIntegrator {

                public:
                    SingleCPUEulerDBIntegrator(SingleCPUKernel *kernel);

                    virtual void execute() override;

                private:
                    SingleCPUKernel *kernel;
                };
            }
        }
    }
}


#endif //READDY_MAIN_SingleCPUEulerDBIntegrator_H
