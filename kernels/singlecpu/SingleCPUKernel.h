//
// Created by clonker on 07.03.16.
//

#ifndef READDY2_MAIN_SINGLECPUKERNEL_H
#define READDY2_MAIN_SINGLECPUKERNEL_H

#include <readdy/plugin/Kernel.h>
#include <boost/dll.hpp>

#define BOOST_DLL_FORCE_ALIAS_INSTANTIATION

namespace readdy {
    namespace kernel {
        namespace singlecpu {
            class SingleCPUKernel : public readdy::plugin::Kernel {
            public:
                SingleCPUKernel();
                ~SingleCPUKernel();
                // move
                SingleCPUKernel(SingleCPUKernel &&rhs);
                SingleCPUKernel& operator=(SingleCPUKernel &&rhs);
                // copy
                SingleCPUKernel(const SingleCPUKernel &rhs);
                SingleCPUKernel& operator=(SingleCPUKernel &rhs);
                // factory method
                static std::shared_ptr<SingleCPUKernel> create();

                virtual std::shared_ptr<readdy::plugin::Program> createProgram(std::string name) const override;
            private:
                struct Impl;
                std::unique_ptr<Impl> pimpl;
            };

            // export factory method as "create_kernel"
            BOOST_DLL_ALIAS(
                    readdy::kernel::singlecpu::SingleCPUKernel::create,
                    createKernel
            );
        }
    }
}

#endif //READDY2_MAIN_SINGLECPUKERNEL_H
