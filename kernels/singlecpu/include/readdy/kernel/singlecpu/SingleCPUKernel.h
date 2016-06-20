//
// Created by clonker on 07.03.16.
//

#ifndef READDY_MAIN_SINGLECPUKERNEL_H
#define READDY_MAIN_SINGLECPUKERNEL_H

#include <readdy/model/RandomProvider.h>
#include <readdy/model/Kernel.h>
#include <boost/dll.hpp>
#include <readdy/kernel/singlecpu/SingleCPUKernelStateModel.h>

#define BOOST_DLL_FORCE_ALIAS_INSTANTIATION

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

                virtual readdy::model::KernelStateModel& getKernelStateModel() const override;
                readdy::kernel::singlecpu::SingleCPUKernelStateModel& getKernelStateModelSingleCPU() const;

                virtual readdy::model::KernelContext& getKernelContext() const override;

                readdy::model::RandomProvider& getRandomProvider() const;

                virtual readdy::model::programs::ProgramFactory &getProgramFactory() const override;

                virtual std::vector<std::string> getAvailablePotentials() const override;

                virtual std::unique_ptr<readdy::model::potentials::Potential> createPotential(std::string &name) const override;

                virtual readdy::model::potentials::PotentialFactory& getPotentialFactory() const override;

            private:
                struct Impl;
                std::unique_ptr<readdy::kernel::singlecpu::SingleCPUKernel::Impl> pimpl;
            };

            BOOST_DLL_ALIAS(
                    readdy::kernel::singlecpu::SingleCPUKernel::name,
                    name
            );
            // export factory method as "create_kernel"
            BOOST_DLL_ALIAS(
                    readdy::kernel::singlecpu::SingleCPUKernel::create,
                    createKernel
            );
        }
    }
}

#endif //READDY_MAIN_SINGLECPUKERNEL_H
