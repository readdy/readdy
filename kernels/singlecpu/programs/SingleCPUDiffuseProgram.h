/**
 * << detailed description >>
 *
 * @file SingleCPUDiffuseProgram.h
 * @brief << brief description >>
 * @author clonker
 * @date 19.04.16
 */

#ifndef READDY_MAIN_SINGLECPUDIFFUSEPROGRAM_H
#define READDY_MAIN_SINGLECPUDIFFUSEPROGRAM_H

#include <readdy/model/Programs.h>
#include "../SingleCPUKernel.h"

namespace readdy {
    namespace kernel {
        namespace singlecpu {

            namespace programs {
                class SingleCPUDiffuseProgram : public readdy::model::DiffuseProgram{

                public:
                    SingleCPUDiffuseProgram(SingleCPUKernel *kernel);

                    virtual ~SingleCPUDiffuseProgram() override = default;

                    virtual void execute() override;

                private:
                    SingleCPUKernel *kernel;
                };
            }
        }
    }
}


#endif //READDY_MAIN_SINGLECPUDIFFUSEPROGRAM_H
