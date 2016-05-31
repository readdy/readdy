//
// Created by clonker on 08.04.16.
//

#ifndef READDY_MAIN_SINGLECPUPROGRAMFACTORY_H
#define READDY_MAIN_SINGLECPUPROGRAMFACTORY_H

#include <readdy/model/ProgramFactory.h>
#include "SingleCPUKernel.h"
#include "SingleCPUKernelStateModel.h"

namespace readdy {
    namespace kernel {
        namespace singlecpu {
            class SingleCPUKernel;
            class SingleCPUProgramFactory : public readdy::model::ProgramFactory {
            public:
                SingleCPUProgramFactory(SingleCPUKernel *kernel);

                virtual std::unique_ptr<readdy::model::Program> createProgram(const std::string& name) const override;

            private:
                SingleCPUKernel *kernel;
            };
        }
    }
}


#endif //READDY_MAIN_SINGLECPUPROGRAMFACTORY_H
