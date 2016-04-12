//
// Created by clonker on 08.04.16.
//

#ifndef READDY2_MAIN_SINGLECPUPROGRAMFACTORY_H
#define READDY2_MAIN_SINGLECPUPROGRAMFACTORY_H

#include <readdy/plugin/ProgramFactory.h>

namespace readdy {
    namespace kernel {
        namespace singlecpu {
            class SingleCPUProgramFactory : readdy::plugin::ProgramFactory<SingleCPUProgramFactory> {
            public:
                static std::shared_ptr<readdy::plugin::Program> createProgram(const std::string name);
            };
        }
    }
}


#endif //READDY2_MAIN_SINGLECPUPROGRAMFACTORY_H
