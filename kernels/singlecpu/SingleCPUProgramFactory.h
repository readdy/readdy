//
// Created by clonker on 08.04.16.
//

#ifndef READDY2_MAIN_SINGLECPUPROGRAMFACTORY_H
#define READDY2_MAIN_SINGLECPUPROGRAMFACTORY_H

#include <readdy/plugin/ProgramFactory.h>
#include "SingleCPUKernel.h"
#include "SingleCPUKernelStateModel.h"

namespace readdy {
    namespace kernel {
        namespace singlecpu {
            class SingleCPUProgramFactory : public readdy::plugin::ProgramFactory {
            public:
                SingleCPUProgramFactory(std::shared_ptr<readdy::model::KernelContext> context, std::shared_ptr<SingleCPUKernelStateModel> model, std::shared_ptr<readdy::utils::RandomProvider> randomProvider);
                virtual std::shared_ptr<readdy::plugin::Program> createProgram(const std::string name) override;

            private:
                std::shared_ptr<readdy::model::KernelContext> context;
                std::shared_ptr<SingleCPUKernelStateModel> model;
                std::shared_ptr<readdy::utils::RandomProvider> randomProvider;
            };
        }
    }
}


#endif //READDY2_MAIN_SINGLECPUPROGRAMFACTORY_H
