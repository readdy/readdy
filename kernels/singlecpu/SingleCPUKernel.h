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
                ~SingleCPUKernel() {
                    BOOST_LOG_TRIVIAL(debug) << "destroying kernel " << getName();
                }
                // factory method
                static std::shared_ptr<SingleCPUKernel> create();
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
