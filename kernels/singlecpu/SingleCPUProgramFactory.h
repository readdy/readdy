//
// Created by clonker on 08.04.16.
//

#ifndef READDY2_MAIN_SINGLECPUPROGRAMFACTORY_H
#define READDY2_MAIN_SINGLECPUPROGRAMFACTORY_H

#include <readdy/plugin/ProgramFactory.h>
#include "SingleCPUKernel.h"

namespace readdy {
    namespace kernel {
        namespace singlecpu {
            class SingleCPUProgramFactory : public readdy::plugin::ProgramFactory {
            public:
                SingleCPUProgramFactory(const SingleCPUKernel &kernel);
                virtual std::shared_ptr<readdy::plugin::Program> createProgram(const std::string name) override;

            private:
                const SingleCPUKernel kernel;
            };
        }
    }
}


#endif //READDY2_MAIN_SINGLECPUPROGRAMFACTORY_H
