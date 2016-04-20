/**
 * << detailed description >>
 *
 * @file SingleCPUDiffuseProgram.h
 * @brief << brief description >>
 * @author clonker
 * @date 19.04.16
 */

#ifndef READDY2_MAIN_SINGLECPUDIFFUSEPROGRAM_H
#define READDY2_MAIN_SINGLECPUDIFFUSEPROGRAM_H

#include <readdy/plugin/Programs.h>
#include "../SingleCPUKernel.h"

namespace readdy {
    namespace kernel {
        namespace singlecpu {
            namespace programs {
                class SingleCPUDiffuseProgram : public readdy::plugin::DiffuseProgram{

                public:
                    SingleCPUDiffuseProgram(std::shared_ptr<SingleCPUKernel> &kernel);

                    virtual ~SingleCPUDiffuseProgram() override = default;

                    virtual void execute() override;

                private:
                    std::shared_ptr<SingleCPUKernel> kernel;
                };
            }
        }
    }
}


#endif //READDY2_MAIN_SINGLECPUDIFFUSEPROGRAM_H
