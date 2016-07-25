//
// Created by clonker on 07.03.16.
//

#ifndef READDY_MAIN_SINGLECPUKERNEL_H
#define READDY_MAIN_SINGLECPUKERNEL_H

#include <readdy/model/RandomProvider.h>
#include <readdy/model/Kernel.h>
#include <readdy/kernel/singlecpu/SingleCPUKernelStateModel.h>

namespace readdy {
    namespace kernel {
        namespace singlecpu {
            class SingleCPUKernel : public readdy::model::Kernel{
            public:

                static const std::string name;

                SingleCPUKernel();
                ~SingleCPUKernel();
                // move
                SingleCPUKernel(SingleCPUKernel &&rhs);
                SingleCPUKernel& operator=(SingleCPUKernel&& rhs);
                // factory method
                static std::unique_ptr<SingleCPUKernel> create();

                virtual SingleCPUKernelStateModel& getKernelStateModel() const override;

                virtual readdy::model::KernelContext& getKernelContext() const override;

                readdy::model::RandomProvider& getRandomProvider() const;

                virtual readdy::model::programs::ProgramFactory &getProgramFactory() const override;

                virtual std::vector<std::string> getAvailablePotentials() const override;

                virtual std::unique_ptr<readdy::model::potentials::Potential> createPotential(std::string &name) const override;

                virtual readdy::model::potentials::PotentialFactory& getPotentialFactory() const override;

                virtual readdy::model::reactions::ReactionFactory &getReactionFactory() const override;

                virtual readdy::model::_internal::ObservableFactory &getObservableFactory() const override;

            private:
                struct Impl;
                std::unique_ptr<readdy::kernel::singlecpu::SingleCPUKernel::Impl> pimpl;
            };
        }
    }
}

#endif //READDY_MAIN_SINGLECPUKERNEL_H
